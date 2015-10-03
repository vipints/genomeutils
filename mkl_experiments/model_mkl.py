#!/usr/bin/env python
"""
mkl version of individual models
"""

from model import ShogunPredictor

class ModelIndividual(object):
    """
    implement light-weight individual using single-task model
    """

    def __init__(self, **kwargs):

        self.task_names = None
        self.kwargs = kwargs
        self.models = {}


    def train(self, data):
        """
        here, we assume data to be dict

        data["C_elegans"]["examples"]
        data["C_elegans"]["labels"]
        """

        self.task_names = data.keys()

        for org in data.keys():
            self.models[org] = ShogunPredictor(**self.kwargs).train(data[org]["examples"], data[org]["labels"])
 

    def predict(self, data):
        """
        predict for each task separately
        """

        return dict((key, self.models[key].predict(data[key]["examples"])) for key in data.keys())

