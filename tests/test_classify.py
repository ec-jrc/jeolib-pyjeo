"""Test suite for module pyjeo.classify."""

import pyjeo as pj
import numpy as np
import unittest


testFile = 'tests/data/s2.tif'
reference = 'tests/data/clc_32632.tif'
model='tests/data/sml.txt'

class BadClassify(unittest.TestCase):
    """Test functions and methods from clasisfy module."""

    def test_classify(self):
        """Test the Symbolic Machine Learning classifier."""

        classDict={}
        classDict['urban']=2
        classDict['agriculture']=12
        classDict['forest']=25
        classDict['water']=41
        classDict['rest']=50
        classFrom=range(0,50)
        classTo=[50]*50
        for i in range(0,50):
            if i>=1 and i<10:
                classTo[i]=classDict['urban']
            elif i>=11 and i<22:
                classTo[i]=classDict['agriculture']
            elif i>=23 and i<25:
                classTo[i]=classDict['forest']
            elif i>=40 and i<45:
                classTo[i]=classDict['water']
            else:
                classTo[i]=classDict['rest']

        jim = pj.Jim(testFile,band2plane=True)
        jim_ref=pj.Jim(reference,dx=jim.properties.getDeltaX(),dy=jim.properties.getDeltaY(),ulx=jim.properties.getUlx(),uly=jim.properties.getUly(),lrx=jim.properties.getLrx(),lry=jim.properties.getLry(),t_srs=jim.properties.getProjection())
        jim_ref.classify.reclass(classes=list(classFrom),reclasses=classTo)

        reflist=pj.JimList([jim_ref])
        jim.classify.trainSML(reflist,output=model,classes=sorted(classDict.values()))
        sml=pj.classify.classify(jim,method='sml',model=model)
        sml.geometry.band2plane()
        sml.np()[:]=np.argmax(sml.np(),axis=0)
        sml.properties.clearNoData()
        sml.classify.reclass(classes=[0,1,2,3,4],reclasses=[2,12,25,41,50])
        stats=sml.stats.getStats('histogram')
        assert stats['histogram'][stats['bin'].index(2)] == 153351.0, \
            'Error in class 0'
        assert stats['histogram'][stats['bin'].index(12)] == 741058.0, \
            'Error in class 1'
        assert stats['histogram'][stats['bin'].index(25)] == 174002.0, \
            'Error in class 2'
        assert stats['histogram'][stats['bin'].index(41)] == 24778.0, \
            'Error in class 3'
        assert stats['histogram'][stats['bin'].index(50)] == 112415.0, \
            'Error in class 4'
        reflist=pj.JimList([jim_ref])


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadClassify)]
    return unittest.TestSuite(suite_list)
