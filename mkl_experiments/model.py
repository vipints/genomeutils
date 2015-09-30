#!/usr/bin/env python
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Written (W) 2010-2013 Christian Widmer
# Copyright (C) 2010-2013 Max-Planck-Society



from shogun.Classifier import SVMLight, LibLinear, MSG_DEBUG, MSG_INFO

from shogun_factory import split_data_promoter, create_promoter_kernel, create_promoter_features 
from shogun_factory import create_labels, create_hashed_features_wdk, create_hashed_promoter_features



class ShogunPredictor(object):
    """
    basic single-task promoter model using string kernels
    """

    def __init__(self, degree=4, shifts=32, kernel_cache=10000, cost=1.0):
        #TODO: clean up degree
        self.degree = degree
        self.degree_wdk = degree
        self.degree_spectrum = degree
        self.shifts = shifts
        self.kernel_cache = kernel_cache
        self.cost = cost
        self.center_offset = 50
        self.center_pos = 1200
        self.epsilon = 10e-2
        self.num_threads = 4


    def train(self, data, labels):

        kernel = create_promoter_kernel(data, self.center_offset, self.center_pos, self.degree_wdk, self.degree_spectrum, self.shifts, kernel_cache=self.kernel_cache)

        print "len(labels) = %i" % (len(labels))
        lab = create_labels(labels)
        self.svm = SVMLight(self.cost, kernel, lab)

        # show debugging output
        self.svm.io.enable_progress()
        self.svm.io.set_loglevel(MSG_DEBUG)

        # optimization settings
        num_threads = self.num_threads
        self.svm.parallel.set_num_threads(num_threads)
        self.svm.set_epsilon(self.epsilon)

        self.svm.train()

        return self


    def predict(self, data):

        feat = create_promoter_features(data, self.center_offset, self.center_pos)
        out = self.svm.apply(feat).get_values()

        return out



