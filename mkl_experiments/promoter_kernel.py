#!/usr/bin/env python
"""
wrapped string kernels for promoter prediction tasks
"""

import uuid
import numpy

from shogun.Kernel import WeightedDegreeStringKernel, CombinedKernel, WeightedCommWordStringKernel, WeightedDegreePositionStringKernel
from shogun.Features import StringCharFeatures, BinaryLabels, DNA, CombinedFeatures, SortWordString, StringWordFeatures
from shogun.Classifier import SVMLight, MSG_DEBUG


def get_spectrum_features(data, order=3, gap=0, reverse=True):
    """
    create feature object used by spectrum kernel
    """

    charfeat = StringCharFeatures(data, DNA)
    feat = StringWordFeatures(charfeat.get_alphabet())
    feat.obtain_from_char(charfeat, order-1, order, gap, reverse)
    preproc = SortWordString()                                            
    preproc.init(feat)
    feat.add_preprocessor(preproc)
    feat.apply_preprocessor()

    return feat


def split_data_promoter(data, center_offset, center_pos):
    """
    divide the sequence in different parts
    """

    center = [seq[(center_pos - center_offset):(center_pos + center_offset)] for seq in data]
    left = [seq[0:center_pos] for seq in data]
    right = [seq[center_pos:] for seq in data]

    return (center, left, right)


def create_features(data, center_offset, center_pos):
    """
    combined features from different kernel perspective 
    """
    
    #print("creating promoter features - center_offset %d, center_pos %d - %s" % (center_offset, center_pos, data))

    (center, left, right) = split_data_promoter(data, center_offset, center_pos)
    # sanity check sequences
    assert len(center) == len(left) == len(right)

    feat_center = StringCharFeatures(DNA)
    feat_center.set_features(center)
    feat_left = get_spectrum_features(left)
    feat_right = get_spectrum_features(right)

    # contruct combined features
    feat = CombinedFeatures()
    feat.append_feature_obj(feat_center)
    feat.append_feature_obj(feat_left)
    feat.append_feature_obj(feat_right)

    return feat


class ShogunPredictor(object):
    """
    basic promoter model using string kernels
    """

    def __init__(self, param):
        self.param = param


    def train(self, data, labels):
        """
        model training 
        """

        # centered WDK/WDK-shift
        if self.param["shifts"] == 0:
            kernel_center = WeightedDegreeStringKernel(self.param["degree"])
        else:
            kernel_center = WeightedDegreePositionStringKernel(10, self.param["degree"])
            shifts_vector = numpy.ones(self.param["center_offset"]*2, dtype=numpy.int32)*self.param["shifts"]
            kernel_center.set_shifts(shifts_vector)

        kernel_center.set_cache_size(self.param["kernel_cache"]/3)

        # border spetrum kernels
        size = self.param["kernel_cache"]/3
        use_sign = False
        kernel_left = WeightedCommWordStringKernel(size, use_sign)
        kernel_right = WeightedCommWordStringKernel(size, use_sign)
        
        # assemble combined kernel
        kernel = CombinedKernel()
        kernel.append_kernel(kernel_center)
        kernel.append_kernel(kernel_left)
        kernel.append_kernel(kernel_right)

        ## building features 
        feat = create_features(data, self.param["center_offset"], self.param["center_pos"])
        
        # init combined kernel
        kernel.init(feat, feat)

        print "len(labels) = %i" % (len(labels))
        lab = BinaryLabels(numpy.double(labels))
        self.svm = SVMLight(self.param["cost"], kernel, lab)

        # show debugging output
        self.svm.io.enable_progress()
        self.svm.io.set_loglevel(MSG_DEBUG)

        # optimization settings
        num_threads = 2
        self.svm.parallel.set_num_threads(num_threads)
        self.svm.set_epsilon(10e-8)

        self.svm.train()

        return self


    def predict(self, data):
        """
        model prediction 
        """
        
        feat = create_features(data, self.param["center_offset"], self.param["center_pos"])
        out = self.svm.apply(feat).get_values()

        return out



