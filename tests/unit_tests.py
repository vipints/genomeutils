import sys 
import logging
import os.path
import unittest

from fetch_remote_data import download_data 


class TestDataProcs(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        sys.stdout.write("unit test")
         
    def testFetchEnsemblFasta(self):

        self.data = {}

         
    def testFetchEnsemblGtf(self):

        self.data = {}


def getTestSuite():
    """
    set up composite test suite
    """
    suite1 = unittest.TestLoader().loadTestsFromTestCase(TestDataProcs)
    return unittest.TestSuite([suite1,suite2])

if __name__ == '__main__':
    unittest.main()
