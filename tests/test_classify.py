"""Test suite for module pyjeo.classify."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
# Copyright (C) 2018-2020 European Union (Joint Research Centre)
#
# This file is part of pyjeo.
#
# pyjeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyjeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyjeo.  If not, see <https://www.gnu.org/licenses/>.

import pyjeo as pj
import numpy as np
import unittest

import os


testFile = 'tests/data/modis_ndvi_2010.tif'
clc = 'tests/data/clc_32632.tif'
trainingfn = 'tests/data/modis_ndvi_training.sqlite'
model = pj._get_random_path()


class BadClassify(unittest.TestCase):
    """Test functions and methods from clasisfy module."""

    @staticmethod
    def test_sml():
        """Test the Symbolic Machine Learning classifier."""
        class_dict = {'urban': 2, 'agriculture': 12, 'forest': 25,
                      'water': 41, 'rest': 50}
        class_from = range(0, 50)
        class_to = [50] * 50
        for i in range(0, 50):
            if 1 <= i < 10:
                class_to[i] = class_dict['urban']
            elif 11 <= i < 22:
                class_to[i] = class_dict['agriculture']
            elif 23 <= i < 25:
                class_to[i] = class_dict['forest']
            elif 40 <= i < 45:
                class_to[i] = class_dict['water']
            else:
                class_to[i] = class_dict['rest']

        jim_ref = pj.Jim(clc, dx=1000, dy=1000)
        jim_ref.classify.reclass(classes=list(class_from), reclasses=class_to)

        bbox = [4246000, 2547000, 4349500, 2441000]
        jim = pj.Jim(testFile, band2plane=True,
                     ulx=bbox[0], uly=bbox[1], lrx=bbox[2], lry=bbox[3])

        jim_ref.geometry.warp(jim.properties.getProjection())
        reflist = pj.JimList([jim_ref])
        jim.classify.trainSML(reflist, output=model,
                              classes=sorted(class_dict.values()))
        sml = pj.classify.sml(jim, model=model)
        sml.geometry.band2plane()
        sml.np()[:] = np.argmax(sml.np(), axis=0)
        sml.properties.clearNoData()
        sml.classify.reclass(classes=[0, 1, 2, 3, 4],
                             reclasses=[2, 12, 25, 41, 50])

        stats = sml.stats.getStats('histogram')
        assert 6605.0 < stats['histogram'][stats['bin'].index(2)] < 6610.0,\
            'Error in class 2'
        assert 23500.0 < stats['histogram'][stats['bin'].index(12)] < 23510.0,\
            'Error in class 12'
        assert 6285.0 < stats['histogram'][stats['bin'].index(25)] < 6295.0,\
            'Error in class 25'
        assert 960.0 < stats['histogram'][stats['bin'].index(41)] < 975.0,\
            'Error in class 41'
        assert 6510.0 < stats['histogram'][stats['bin'].index(50)] < 6520.0,\
            'Error in class 50'
        os.remove(model)


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadClassify)]
    return unittest.TestSuite(suite_list)
