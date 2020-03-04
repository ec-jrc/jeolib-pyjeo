"""Test suite for module pyjeo.classify."""

import pyjeo as pj
import numpy as np
import unittest

import os


testFile = 'tests/data/modis_ndvi_2010.tif'
reference = 'tests/data/clc_32632.tif'
model = pj._get_random_path()


class BadClassify(unittest.TestCase):
    """Test functions and methods from clasisfy module."""

    def test_classify(self):
        """Test the Symbolic Machine Learning classifier."""
        classDict = {}
        classDict['urban'] = 2
        classDict['agriculture'] = 12
        classDict['forest'] = 25
        classDict['water'] = 41
        classDict['rest'] = 50
        classFrom = range(0, 50)
        classTo = [50] * 50
        for i in range(0, 50):
            if 1 <= i < 10:
                classTo[i] = classDict['urban']
            elif 11 <= i < 22:
                classTo[i] = classDict['agriculture']
            elif 23 <= i < 25:
                classTo[i] = classDict['forest']
            elif 40 <= i < 45:
                classTo[i] = classDict['water']
            else:
                classTo[i] = classDict['rest']

        jim_ref = pj.Jim(reference, dx=1000, dy=1000)
        jim_ref.classify.reclass(classes=list(classFrom), reclasses=classTo)

        bbox = [4246000, 2547000, 4349500, 2441000]
        jim = pj.Jim(testFile, band2plane=True,
                     ulx=bbox[0], uly=bbox[1], lrx=bbox[2], lry=bbox[3])

        jim_ref.geometry.warp(jim.properties.getProjection())
        reflist = pj.JimList([jim_ref])
        jim.classify.trainSML(reflist, output=model,
                              classes=sorted(classDict.values()))
        sml = pj.classify.classify(jim, method='sml', model=model)
        sml.geometry.band2plane()
        sml.np()[:] = np.argmax(sml.np(), axis=0)
        sml.properties.clearNoData()
        sml.classify.reclass(classes=[0, 1, 2, 3, 4],
                             reclasses=[2, 12, 25, 41, 50])

        stats = sml.stats.getStats('histogram')
        assert stats['histogram'][stats['bin'].index(2)] == 6608.0, \
            'Error in class 2'
        assert stats['histogram'][stats['bin'].index(12)] == 23507.0, \
            'Error in class 12'
        assert stats['histogram'][stats['bin'].index(25)] == 6289.0, \
            'Error in class 25'
        assert stats['histogram'][stats['bin'].index(41)] == 966.0, \
            'Error in class 41'
        assert stats['histogram'][stats['bin'].index(50)] == 6514.0, \
            'Error in class 50'

        os.remove(model)


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadClassify)]
    return unittest.TestSuite(suite_list)
