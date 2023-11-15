"""Test suite for module pyjeo.geometry."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
# Copyright (C) 2018-2023 European Union (Joint Research Centre)
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
import unittest
import warnings
import numpy as np
# from datetime import time, timedelta, datetime
from datetime import datetime
import os

tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']

rasterfn = 'tests/data/modis_ndvi_2010.tif'
vectorfn = 'tests/data/modis_ndvi_training.sqlite'
nutsfn = 'tests/data/nuts_italy.sqlite'
clc = 'tests/data/clc_32632.tif'
outputfn = pj._get_random_path()

class BadGeometry(unittest.TestCase):
    """Test functions and methods from geometry module."""

    @staticmethod
    def test_covers():
        """Test the covers method."""
        jim= pj.Jim(rasterfn, band=0)
        bbox = jim.properties.getBBox()
        ulx = bbox [0]
        uly = bbox [1]
        lrx = bbox [2]
        lry = bbox [3]
        bbox = [ulx, uly, lrx, lry]
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry,
                                coverType = 'ALL_COVERED') ==
            jim.geometry.covers(bbox = bbox,
                                coverType = 'ALL_COVERED')), \
            'Error in method geometry.covers() bbox'
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry,
                                coverType = 'ALL_TOUCHED') ==
            jim.geometry.covers(bbox = bbox, coverType = 'ALL_TOUCHED')), \
            'Error in method geometry.covers() bbox'
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry) ==
                jim.geometry.covers(bbox = bbox)), \
                'Error in method geometry.covers() bbox'
        assert(pj.geometry.covers(jim, ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry,
                                coverType = 'ALL_COVERED') ==
            pj.geometry.covers(jim, bbox = bbox,
                                coverType = 'ALL_COVERED')), \
            'Error in function geometry.covers() bbox'
        assert(pj.geometry.covers(jim, ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry,
                                coverType = 'ALL_TOUCHED') ==
            pj.geometry.covers(jim, bbox = bbox,
                                coverType = 'ALL_TOUCHED')), \
            'Error in function geometry.covers() bbox'
        assert(pj.geometry.covers(jim, ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry) ==
            pj.geometry.covers(jim, bbox = bbox)), \
            'Error in function geometry.covers() bbox'
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry,
                                coverType = 'ALL_COVERED') ==
            pj.geometry.covers(jim, bbox = bbox,
                                coverType = 'ALL_COVERED')), \
            'Error in function geometry.covers() function not equal to method'
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry,
                                coverType = 'ALL_TOUCHED') ==
            pj.geometry.covers(jim, bbox = bbox,
                                coverType = 'ALL_TOUCHED')), \
            'Error in function geometry.covers() function not equal to method'
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                    lrx = lrx, lry = lry) ==
                pj.geometry.covers(jim, bbox = bbox)), \
                'Error in function geometry.covers() function not equal to method'
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry,
                                coverType = 'ALL_COVERED') == False), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry,
                                coverType = 'ALL_TOUCHED') == True), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = ulx, uly = uly,
                                lrx = lrx, lry = lry) == True), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = ulx + 10, uly = uly - 10,
                                lrx = lrx - 10, lry = lry + 10,
                                coverType = 'ALL_COVERED') == True), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = ulx + 10, uly = uly - 10,
                                lrx = lrx - 10, lry = lry + 10,
                                coverType = 'ALL_TOUCHED') == True), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = ulx + 10, uly = uly - 10,
                                lrx = lrx - 10, lry = lry + 10) == True), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = ulx + 10, uly = uly - 10,
                                lrx = lrx + 10, lry = lry - 10,
                                coverType = 'ALL_COVERED') == False), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = ulx + 10, uly = uly - 10,
                                lrx = lrx - 10, lry = lry + 10,
                                coverType = 'ALL_TOUCHED') == True), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = ulx + 10, uly = uly - 10,
                                lrx = lrx - 10, lry = lry + 10) == True), \
                'Error in method geometry.covers() must be True'
        assert(jim.geometry.covers(ulx = 2*ulx, uly = uly,
                                lrx = 2*lrx, lry = lry,
                                coverType = 'ALL_COVERED') == False), \
                                'Error in method geometry.covers() must be False'
        assert(jim.geometry.covers(ulx = 2*ulx, uly = uly,
                                lrx = 2*lrx, lry = lry,
                                coverType = 'ALL_TOUCHED') == False), \
                                'Error in method geometry.covers() must be False'
        assert(jim.geometry.covers(ulx = 2*ulx, uly = uly,
                                lrx = 2*lrx, lry = lry) == False), \
                                'Error in method geometry.covers() must be False'
        assert(jim.geometry.covers(ulx = 2*ulx, uly = uly,
                                lrx = 2*lrx, lry = lry,
                                coverType = 'ALL_COVERED') == False), \
                                'Error in method geometry.covers() must be False'

    @staticmethod
    def test_band2plane():
        """Test the band2plane method."""
        jim3d = pj.Jim(rasterfn, band=[0, 1], band2plane=True)
        jim2d = pj.Jim(rasterfn, band=[0, 1], band2plane=False)
        jim2d.geometry.band2plane()
        assert jim2d.properties.isEqual(jim3d), \
            'Error in geometry.band2plane() ' \
            '(read as 3d is not equal to convert to 3d)'
        jim2d.geometry.plane2band()
        assert pj.geometry.plane2band(jim3d).properties.isEqual(jim2d), \
            'Error in geometry.plane2band() ' \
            '(function is not equal to method)'
        assert pj.geometry.band2plane(jim2d).properties.isEqual(jim3d), \
            'Error in geometry.band2plane() ' \
            '(function is not equal to method)'

    @staticmethod
    def test_band2plane_dimension():
        jim3d = pj.Jim(rasterfn, band=[0, 1], band2plane=True)
        dates = [
            datetime.strptime('2019-01-01','%Y-%m-%d'),
            datetime.strptime('2019-01-05','%Y-%m-%d')]
        jim3d.properties.setDimension({'plane' : dates,
                                        'band' : ['B2']})
        jim2d = pj.Jim(rasterfn, band=[0, 1], band2plane=False)
        jim2d.geometry.band2plane()
        jim2d.properties.setDimension({'plane' : dates,
                                        'band':['B3']})
        assert jim2d.properties.isEqual(jim3d), \
            'Error in geometry.band2plane() ' \
            '(read as 3d is not equal to convert to 3d)'

        jim2d.geometry.plane2band()
        assert jim2d.properties.getDimension(
            'band') == jim3d.properties.getDimension('plane'), \
            'Error in geometry.band2plane() dimension band / temporal'

        assert pj.geometry.plane2band(jim3d).properties.getDimension(
            'band') == jim2d.properties.getDimension('band'), \
            'Error in function pj.geometry.plane2band()'

        assert pj.geometry.band2plane(jim2d) == jim3d, \
            'Error in geometry band2plane'

        assert pj.geometry.band2plane(jim2d).properties.getDimension(
            'plane') == jim3d.properties.getDimension('plane'), \
            'Error in geometry band2plane, dimension'


    @staticmethod
    def test_stack():
        """Test the stackBand and stackPlane methods."""
        jim = pj.Jim(rasterfn, band=[0, 1, 2])

        # Test stackBand()
        jim0to1 = pj.geometry.cropBand(jim, [0, 1])
        jim2 = pj.geometry.cropBand(jim, 2)

        jim_stacked = pj.geometry.stackBand(jim0to1, jim2)
        assert jim_stacked.properties.nrOfBand() == \
               jim0to1.properties.nrOfBand() + jim2.properties.nrOfBand(), \
            'Error in geometry.stackBand() ' \
            '(number of bands after stack not equal to the sum of number of ' \
            'bands of parameters)'

        jim0to1.geometry.stackBand(jim2)

        assert jim_stacked.properties.isEqual(jim0to1), \
            'Inconsistency in geometry.stackBand() ' \
            '(method returns different result than function)'
        assert jim_stacked.properties.isEqual(jim), \
            'Error in geometry.stackBand() ' \
            '(stacked bands are not equal to an object from which they were ' \
            'cropped)'

        # Test stackBand() with band specified
        jim0to1 = pj.geometry.cropBand(jim, [0, 1])

        jim_stacked = pj.geometry.stackBand(jim0to1, jim, band=2)
        assert jim_stacked.properties.nrOfBand() == \
               jim.properties.nrOfBand(), \
            'Error in geometry.stackBand(band) ' \
            '(number of bands after stack not equal to the expected one)'

        jim0to1.geometry.stackBand(jim, band=2)

        assert jim_stacked.properties.isEqual(jim0to1), \
            'Inconsistency in geometry.stackBand(band) ' \
            '(method returns different result than function)'
        assert jim_stacked.properties.isEqual(jim), \
            'Error in geometry.stackBand(band) ' \
            '(stacked bands are not equal to an object from which' \
            'they were cropped)'

        # Test the stackPlane() method
        jim0 = pj.Jim(rasterfn, band=0)
        jim1 = pj.Jim(rasterfn, band=1)
        jim2 = pj.Jim(rasterfn, band=2)
        jimlist = pj.JimList([jim0, jim1, jim2])
        jimliststack = pj.geometry.stackPlane(jimlist)
        jimstack = pj.geometry.stackPlane(pj.geometry.stackPlane(jim0, jim1),
                                          jim2)
        assert jimliststack.properties.isEqual(jimstack), \
            'Error in geometry.stackPlane() ' \
            '(jimliststack not equal to jimstack)'
        jim3 = pj.JimList([jim0, jim1, jim2]).geometry.stackPlane()
        jim3.geometry.cropPlane([0, 1])
        assert pj.geometry.cropPlane(jim3, 0).properties.isEqual(jim0), \
            'Error jim3 not equal to jim0'
        jim4 = pj.JimList([jim0, jim1, jim2]).geometry.stackPlane()
        jim4.geometry.cropPlane([1, 2])
        assert pj.geometry.cropPlane(jim4, 0).properties.isEqual(jim1), \
            'Error jim3 not equal to jim1'
        jim = pj.JimList([jim3, jim4]).geometry.stackPlane()
        jim.geometry.cropPlane([1, 2])
        jim.geometry.cropPlane(0)
        assert jim.properties.isEqual(jim1), \
            'Error jim not equal to jim1'

        # Test stackPlane() with only one Jim as a parameter
        stacked = pj.geometry.stackPlane(jim)

        assert stacked.properties.isEqual(jim), \
            'Error in geometry.stackPlane(Jim) ' \
            '(returned object not equal to the only Jim used as a parameter)'

        jim.geometry.stackPlane()

        assert stacked.properties.isEqual(jim), \
            'Inconsistency in geometry.stackPlane(Jim) ' \
            '(method returns different result than function)'

        # Test the stackPlane() with self

        jim1 = pj.Jim(rasterfn)
        jim2 = pj.Jim(jim1)
        jim1.geometry.stackPlane(pj.Jim(jim1))
        jim2.geometry.stackPlane(jim2)
        assert(jim1.properties.isEqual(jim2)), \
            'Inconsistency in geometry.stackPlane(self)'

        # Test the band2plane() method
        jimband = pj.Jim(rasterfn, band=[0, 1, 2])
        jimplane = pj.Jim(rasterfn, band=[0, 1, 2], band2plane=True)
        jimband.geometry.band2plane()
        assert jimband.properties.isEqual(jimplane), \
            'Error in geometry.band2plane() ' \
            '(jimband not equal to jimsplane)'

        # Test wrong call of stackBand with the parameter band exceeding
        # the nrOfBand in Jims
        try:
            _ = pj.geometry.stackBand(jim0, jim1,
                                      band=jim0.properties.nrOfBand() + 1)
            raised = False
        except pj.exceptions.JimPlanesError:
            raised = True

            assert raised, \
                'Error in catching a call of geometry.stackBand(jim, jim, band) ' \
                'function where planes do not match'

        # Test wrong call of stackBand with the parameter band exceeding
        # the nrOfBand in Jims
        jim1.geometry.cropPlane(0)
        try:
            _ = pj.geometry.stackBand(jim0, jim1,
                                    band = jim1.properties.nrOfBand())
            raised = False
        except pj.exceptions.JimBandsError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackBand(jim, jim, band) ' \
            'function where the band argument exceeds the allowed number of ' \
            'bands of Jim'

        # Test wrong call of stackBand with wrong object type
        try:
            _ = pj.geometry.stackBand(1, jim1)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackBand(jim, jim) ' \
            'function where the first argument is not a Jim or a JimList'

        # Test wrong call of stackPlane with wrong object type
        try:
            _ = pj.geometry.stackPlane(1, jim1)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackPlane(jim, jim) ' \
            'function where the first argument is not a Jim or a JimList'


    @staticmethod
    def test_stack_dimension():
        """Test the stackBand and stackPlane methods."""
        jim = pj.Jim(rasterfn, band=[0, 1, 2])

        dates = [
            datetime.strptime('2019-01-01','%Y-%m-%d'),
            datetime.strptime('2019-02-01','%Y-%m-%d'),
            datetime.strptime('2019-03-01','%Y-%m-%d')]
        jim.properties.setDimension(dates, 'band')

        # Test stackBand() with empty Jim
        jimempty = pj.Jim()
        jim0to1 = pj.geometry.cropBand(jim, [0, 1])
        jimempty.geometry.stackBand(jim0to1)
        assert jimempty.properties.getDimension('band') == \
            dates[0:2], \
            'Error in geometry.stackBand() with empty Jim dimension band'
        assert jimempty.properties.isEqual(jim0to1), \
            'Error in geometry.stackBand() with empty Jim'

        jimempty = pj.Jim()
        jim0to1 = pj.geometry.cropBand(jim, [0, 1])
        jimempty = pj.geometry.stackBand(jimempty, jim0to1)
        assert jimempty.properties.getDimension('band') == \
            dates[0:2], \
            'Error in geometry.stackBand() with empty Jim dimension band'
        assert jimempty.properties.isEqual(jim0to1), \
            'Error in geometry.stackBand() with empty Jim and Jim'

        jimempty = pj.Jim()
        jim0to1 = pj.geometry.cropBand(jim, [0, 1])
        jimempty = pj.geometry.stackBand(jim0to1, jimempty)
        assert jimempty.properties.getDimension('band') == \
            dates[0:2], \
            'Error in geometry.stackBand() with empty Jim dimension band'
        assert jimempty.properties.isEqual(jim0to1), \
            'Error in geometry.stackBand() with Jim and empty Jim'

        jimempty = pj.Jim()
        jimempty.geometry.stackBand(jim0to1, 0)
        assert jimempty.properties.getDimension('band') == \
            dates[0:1], \
            'Error in geometry.stackBand() with empty Jim selected band'
        assert jimempty.properties.isEqual(jim0to1), \
            'Error in geometry.stackBand() with empty Jim'

        jimempty = pj.Jim()
        jimempty.geometry.stackBand(
            jim0to1, datetime.strptime('2019-01-01','%Y-%m-%d'),)
        assert jimempty.properties.getDimension('band') == \
            dates[0:1], \
            'Error in geometry.stackBand() with empty Jim selected band'
        assert jimempty.properties.isEqual(jim0to1), \
            'Error in geometry.stackBand() with empty Jim'

        # Test stackBand()
        jim0to1 = pj.geometry.cropBand(jim, [0, 1])
        assert jim0to1.properties.getDimension('band') == \
            dates[0:2], \
            'Error in geometry.stackBand() dimension band'

        assert jim0to1 == pj.geometry.cropBand(jim, dates[0:2]), \
            'Error in geometry.stackBand(): cropBand int != datetime'
        assert pj.geometry.cropBand(jim, dates[0:2]).properties.\
            getDimension('band') == dates[0:2], \
            'Error in geometry.stackBand() dimension band'

        jim2 = pj.geometry.cropBand(jim, dates[2])
        assert jim2.properties.getDimension('band') == \
            dates[2:3], \
            'Error in geometry.stackBand() dimension band'

        jim_stacked = pj.geometry.stackBand(jim0to1, jim2)
        assert jim_stacked.properties.nrOfBand() == \
               jim0to1.properties.nrOfBand() + jim2.properties.nrOfBand(), \
            'Error in geometry.stackBand() ' \
            '(number of bands after stack not equal to the sum of number of ' \
            'bands of parameters)'

        assert jim_stacked.properties.getDimension('band') == dates, \
            'Error in geometry.stackBand() dimensions'

        jim0to1.geometry.stackBand(jim2)

        assert jim_stacked.properties.isEqual(jim0to1), \
            'Inconsistency in geometry.stackBand() ' \
            '(method returns different result than function)'

        assert jim_stacked.properties.getDimension('band') == \
            jim0to1.properties.getDimension('band'), \
            'Error in geometry.stackBand() dimensions'

        assert jim_stacked.properties.isEqual(jim), \
            'Error in geometry.stackBand() ' \
            '(stacked bands are not equal to an object from which they \
            were cropped)'

        # Test stackBand() with band specified
        jim0to1 = pj.geometry.cropBand(jim, dates[0:2])

        jim_stacked = pj.geometry.stackBand(jim0to1, jim, band=2)
        assert jim_stacked.properties.nrOfBand() == \
               jim.properties.nrOfBand(), \
            'Error in geometry.stackBand(band) ' \
            '(number of bands after stack not equal to the expected one)'

        jim0to1.geometry.stackBand(jim, band=dates[2])

        assert jim_stacked.properties.isEqual(jim0to1), \
            'Inconsistency in geometry.stackBand(band) ' \
            '(method returns different result than function)'
        assert jim_stacked.properties.isEqual(jim), \
            'Error in geometry.stackBand(band) ' \
            '(stacked bands are not equal to an object from which they were ' \
            'cropped)'

        assert jim.properties.getDimension('band') == jim_stacked.\
            properties.getDimension('band'), \
            'Error in geometry.stackBand() dimensions with band specified'

        # Test the stackPlane() with empty Jim

        jim1 = pj.Jim(rasterfn)
        jim1.properties.setDimension(dates[1], 'plane')
        jim2 = pj.Jim()
        jim3 = pj.geometry.stackPlane(jim1, jim2)
        assert jim3.properties.getDimension() == \
            jim1.properties.getDimension(), \
            'Inconsistency in dimension geometry.stackPlane(jim, empty)'
        assert(jim3.properties.isEqual(jim1)), \
            'Inconsistency in geometry.stackPlane(jim, empty)'

        jim1 = pj.Jim(rasterfn)
        jim1.properties.setDimension(dates[1], 'plane')
        jim2 = pj.Jim()
        jim3 = pj.geometry.stackPlane(jim2, jim1)
        assert jim3.properties.getDimension() == \
            jim1.properties.getDimension(), \
            'Inconsistency in dimension geometry.stackPlane(empty, jim)'
        assert(jim3.properties.isEqual(jim1)), \
            'Inconsistency in geometry.stackPlane(empty, jim)'

        jim1 = pj.Jim(rasterfn)
        jim1.properties.setDimension(dates[1], 'plane')
        jim2 = pj.Jim()
        jim3 = pj.geometry.stackPlane(jim1, [jim2])
        assert jim3.properties.getDimension() == \
            jim1.properties.getDimension(), \
            'Inconsistency in dimension geometry.stackPlane(jim, [empty])'
        assert(jim3.properties.isEqual(jim1)), \
            'Inconsistency in geometry.stackPlane(empty, jim)'

        jim1 = pj.Jim(rasterfn)
        jim1.properties.setDimension(dates[1], 'plane')
        jim2 = pj.Jim()
        jim3 = pj.geometry.stackPlane(jim2, [jim1])
        assert jim3.properties.getDimension() == \
            jim1.properties.getDimension(), \
            'Inconsistency in dimension geometry.stackPlane(empty, [jim])'
        assert(jim3.properties.isEqual(jim1)), \
            'Inconsistency in geometry.stackPlane(empty, jim)'

        jim1 = pj.Jim(rasterfn)
        jim1.properties.setDimension(dates[1], 'plane')
        jim2 = pj.Jim()
        jim3 = pj.geometry.stackPlane(jim1, [jim2, pj.Jim()])
        assert jim3.properties.getDimension() == \
            jim1.properties.getDimension(), \
            'Inconsistency in dimension geometry.stackPlane(jim, [empty])'
        assert(jim3.properties.isEqual(jim1)), \
            'Inconsistency in geometry.stackPlane(empty, jim)'

        jim1 = pj.Jim(rasterfn)
        jim1.properties.setDimension(dates[1], 'plane')
        jim2 = pj.Jim()
        jim3 = pj.geometry.stackPlane(jim2, [pj.Jim(), jim1, pj.Jim()])
        assert jim3.properties.getDimension() == \
            jim1.properties.getDimension(), \
            'Inconsistency in dimension geometry.stackPlane(empty, [jim])'
        assert(jim3.properties.isEqual(jim1)), \
            'Inconsistency in geometry.stackPlane(empty, jim)'

        jim1 = pj.Jim(rasterfn)
        jim1.properties.setDimension(dates[1], 'plane')
        jim2 = pj.Jim()
        jim1.geometry.stackPlane(jim2)
        assert jim1.properties.getDimension() == \
            jim1.properties.getDimension(), \
            'Inconsistency in dimension geometry.stackPlane(empty, [jim])'
        assert(jim1.properties.isEqual(jim1)), \
            'Inconsistency in geometry.stackPlane(empty, jim)'

        jim1 = pj.Jim(rasterfn)
        jim1.properties.setDimension(dates[1], 'plane')
        jim2 = pj.Jim()
        jim2.geometry.stackPlane(jim1)
        assert jim2.properties.getDimension() == \
            jim1.properties.getDimension(), \
            'Inconsistency in dimension geometry.stackPlane(empty, [jim])'
        assert(jim2.properties.isEqual(jim1)), \
            'Inconsistency in geometry.stackPlane(empty, jim)'

        # Test the stackPlane() with empty JimList
        jimlist = pj.JimList([])
        jimliststack = pj.geometry.stackPlane(jimlist)
        assert not jimliststack, \
            'Error in pj.geometry.stackPlane(jimlist) function empty list'
        assert jimlist.geometry.stackPlane().properties.isEqual(pj.Jim()), \
            'Error in pj.geometry.stackPlane(jimlist) method empty list'

        # Test the stackPlane() method

        jim0 = pj.Jim(rasterfn, band=0)
        jim1 = pj.Jim(rasterfn, band=1)
        jim2 = pj.Jim(rasterfn, band=2)
        jim0.properties.setDimension(dates[0], 'plane')
        jim1.properties.setDimension(dates[1], 'plane')
        jim2.properties.setDimension(dates[2], 'plane')

        jimlist = pj.JimList([jim0, jim1, jim2])
        jimliststack = pj.geometry.stackPlane(jimlist)

        assert jimliststack.properties.getDimension('plane') == dates, \
            'Error in pj.geometry.stackPlane(jimlist) dimensions'

        jimstack = pj.geometry.stackPlane(pj.geometry.stackPlane(jim0, jim1),
                                          jim2)
        assert jimliststack.properties.isEqual(jimstack), \
            'Error in geometry.stackPlane() ' \
            '(jimliststack not equal to jimstack)'

        assert jimstack.properties.getDimension('plane') == dates, \
            'Error in pj.geometry.stackPlane dimensions'

        jim3 = pj.JimList([jim0, jim1, jim2]).geometry.stackPlane()
        jim3.geometry.cropPlane(dates[0:2])

        assert jim3.properties.getDimension('plane') == dates[0:2], \
            'Error in pj.geometry.cropPlane dimensions'

        assert pj.geometry.cropPlane(jim3, 0).properties.isEqual(jim0), \
            'Error jim3 not equal to jim0'

        assert pj.geometry.cropPlane(jim3, 0).properties.getDimension(
            'plane') == dates[0:1], \
            'Error in pj.geometry.cropPlane dimensions'

        jim4 = pj.JimList([jim0, jim1, jim2]).geometry.stackPlane()
        jim4.geometry.cropPlane(dates[1:3])
        assert pj.geometry.cropPlane(jim4, dates[1]).properties.isEqual(jim1), \
            'Error cropped jim4 not equal to jim1'

        assert jim4.properties.getDimension('plane') == dates[1:3], \
            'Error in pj.geometry.cropPlane dimensions'

        jim = pj.JimList([jim3, jim4]).geometry.stackPlane()
        jim.geometry.cropPlane(dates[1:3])
        jim.geometry.cropPlane(dates[1])
        assert jim.properties.isEqual(jim1), \
            'Error jim not equal to jim1'
        assert jim.properties.getDimension('plane') == dates[1:2], \
            'Error in pj.geometry.cropPlane dimensions'

        # Test stackPlane() with only one Jim as a parameter
        stacked = pj.geometry.stackPlane(jim)

        assert stacked.properties.isEqual(jim), \
            'Error in geometry.stackPlane(Jim) ' \
            '(returned object not equal to the only Jim used as a parameter)'

        jim.geometry.stackPlane()

        assert stacked.properties.isEqual(jim), \
            'Inconsistency in geometry.stackPlane(Jim) ' \
            '(method returns different result than function)'

        # Test the stackPlane() with self

        jim1 = pj.Jim(rasterfn)
        jim2 = pj.Jim(jim1)
        jim1.geometry.stackPlane(pj.Jim(jim1))
        jim2.geometry.stackPlane(jim2)
        assert(jim1.properties.isEqual(jim2)), \
            'Inconsistency in geometry.stackPlane(self)'

        # Test the band2plane() method
        jimband = pj.Jim(rasterfn, band=[0, 1, 2])
        jimplane = pj.Jim(rasterfn, band=[0, 1, 2], band2plane=True)
        jimband.geometry.band2plane()
        assert jimband.properties.isEqual(jimplane), \
            'Error in geometry.band2plane() ' \
            '(jimband not equal to jimsplane)'

        # Test wrong call of stackBand with number of planes not matching
        # the nrOfBand in Jims
        try:
            _ = pj.geometry.stackBand(jim0, jim1,
                                      band = jim0.properties.nrOfBand() - 1)
            raised = False
        except pj.exceptions.JimPlanesError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackBand(jim, jim, band) ' \
            'function where planes do not match'
        # Test wrong call of stackBand with the parameter band exceeding
        # the nrOfBand in Jims
        jim1.geometry.cropPlane(0)

        try:
            _ = pj.geometry.stackBand(jim0, jim1,
                                    band = jim1.properties.nrOfBand())
            raised = False
        except pj.exceptions.JimBandsError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackBand(jim, jim, band) ' \
            'function where the band argument exceeds the allowed number of ' \
            'bands of Jim'

        # Test wrong call of stackBand with wrong object type
        try:
            _ = pj.geometry.stackBand(1, jim1)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackBand(jim, jim) ' \
            'function where the first argument is not a Jim or a JimList'

        # Test wrong call of stackPlane with wrong object type
        try:
            _ = pj.geometry.stackPlane(1, jim1)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackPlane(jim, jim) ' \
            'function where the first argument is not a Jim or a JimList'
    @staticmethod
    def test_warp():
        """Test the warp method."""
        jim0 = pj.Jim(rasterfn, band=[0, 1])
        jim1 = pj.Jim(rasterfn, band=[2, 3])
        jim0.geometry.stackPlane(jim1)
        jim0_2 = pj.Jim(jim0)

        orig_proj = jim0.properties.getProjection()

        jim_warped = pj.geometry.warp(jim0, 'epsg:4326')
        jim0_2.geometry.warp('epsg:4326')

        assert jim_warped.properties.isEqual(jim0_2), \
            'Inconsistency in geometry.warp() ' \
            '(method returns different result than function)'
        assert jim_warped.properties.getProjection() != orig_proj, \
            'Error in geometry.warp() ' \
            '(EPSG not changed)'
        assert jim_warped.properties.getProjection()[-7:-3] == '4326', \
            'Error in geometry.warp() ' \
            '(EPSG not changed to the right one)'

        jim_warped.geometry.warp('epsg:3035')

        assert jim_warped.properties.getProjection()[-7:-3] == '3035', \
            'Error in geometry.warp() ' \
            '(EPSG not changed to the right one when called for the second ' \
            'time)'
        assert jim_warped.properties.nrOfBand() == 2, \
            'Error in geometry.warp() ' \
            '(number of bands changed)'
        assert jim_warped.properties.nrOfPlane() == 2, \
            'Error in geometry.warp() ' \
            '(number of planes changed)'

        jim_warped.geometry.crop(ulx=jim0.properties.getUlx(),
                                 uly=jim0.properties.getUly(),
                                 lrx=jim0.properties.getLrx(),
                                 lry=jim0.properties.getLry(),
                                 dx=jim0.properties.getDeltaX(),
                                 dy=jim0.properties.getDeltaY())

        assert jim0.properties.nrOfCol() == jim_warped.properties.nrOfCol(), \
            'Error in geometry.warp() ' \
            '(after warping to a different projection and back and ' \
            'after cropping to the original extent, nrOfCol changed)'
        assert jim0.properties.nrOfRow() == jim_warped.properties.nrOfRow(), \
            'Error in geometry.warp() ' \
            '(after warping to a different projection and back and ' \
            'after cropping to the original extent, nrOfRow changed)'

        #Test the warp method without re-projection for resampling

        jim0 = pj.Jim(rasterfn, band=[0, 1])

        dates = [
            datetime.strptime('2019-01-01','%Y-%m-%d'),
            datetime.strptime('2019-01-05','%Y-%m-%d')]
        bands = ['B0', 'B1']

        jim0.geometry.stackPlane(pj.Jim(rasterfn, band=[2, 3]))
        jim0.properties.setDimension({'plane' : dates,
                                      'band' : bands})

        bbox = jim0.properties.getBBox()
        dx = jim0.properties.getDeltaX()
        dy = jim0.properties.getDeltaY()
        ulx, uly, lrx, lry = bbox
        projection = jim0.properties.getProjection()
        bbox_extended = [ulx - dx,
                         uly + dy,
                         lrx + dx,
                         lry - dy]

        jim_near = pj.geometry.warp(jim0, bbox = bbox_extended,
                                    resample = 'near',
                                    nodata = 255)
        assert jim_near.np(0)[0,0,0] == 255, \
            'Error in geometry.warp(): extend bounds 0 0 0 with nodata 255 band 0'
        assert jim_near.np(1)[0,0,0] == 255, \
            'Error in geometry.warp(): extend bounds 0 0 0 with nodata 255 band 1'
        assert jim_near.np(0)[-1,-1,-1] == 255, \
            'Error in geometry.warp(): extend bounds -1 -1 -1 with nodata 255 band 0'
        assert jim_near.np(1)[-1,-1,-1] == 255, \
            'Error in geometry.warp(): extend bounds -1 -1 -1 with nodata 255 band 1'
        assert jim_near.properties.getBBox() == bbox_extended, \
            'Error in geometry.warp(): bbox_extended'
        assert jim_near.properties.getProjection() == projection, \
            'Error in geometry.warp(): projection'

        jim_near = pj.geometry.warp(jim0, bbox = bbox,
                                    resample = 'near',
                                    nodata = 255)
        assert jim_near.properties.isEqual(jim0), \
            'Error in geometry.warp(): recover original from extended jim'
        assert jim_near.properties.getBBox() == bbox, \
            'Error in geometry.warp(): bbox'
        assert jim_near.properties.getProjection() == projection, \
            'Error in geometry.warp(): projection'

        #Test the warp method without re-projection for resampling

        projection = jim0.properties.getProjection()

        jim_cubic = pj.geometry.warp(jim0, dx = 100, dy = 100,
                                     resample = 'cubic')

        jim_near = pj.geometry.warp(jim0, dx = 100, dy = 100,
                                    resample = 'near')
        assert not jim_cubic.properties.isEqual(jim_near), \
            'Error in geometry.warp(): resampling method'

        jim0.geometry.warp(dx = 100, dy = 100, resample = 'cubic')

        assert jim_cubic.properties.isEqual(jim0), \
            'Error in geometry.warp() ' \
            'function not equal to method'

        assert jim0.properties.nrOfPlane() == len(dates), \
            'Error in geometry.warp(): len plane'
        assert jim0.properties.nrOfBand() == len(bands), \
            'Error in geometry.warp(): len band'
        assert jim0.properties.getDimension('plane') == dates, \
            'Error in geometry.warp(): dimension plane'
        assert jim0.properties.getDimension('band') == bands, \
            'Error in geometry.warp(): dimension band'
        assert jim0.properties.getDeltaX() == 100, \
            'Error in geometry.warp(): dx is not 100'
        assert jim0.properties.getDeltaY() == 100, \
            'Error in geometry.warp(): dy is not 100'
        assert jim0.properties.getBBox() == bbox, \
            'Error in geometry.warp(): bbox'
        assert jim0.properties.getProjection() == projection, \
            'Error in geometry.warp(): projection'

    @staticmethod
    def test_extract_loop():
        """Test the extract method looping over bands using join function."""
        jim0 = pj.Jim(rasterfn)
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()

        nband = jim0.properties.nrOfBand()

        for band in range(0, nband):
            jim = pj.Jim(pj.geometry.cropBand(jim0, band))
            bandname = ['B' + str(band)]
            if not band:
                v = pj.geometry.extract(sampleid, jim, rule='mean',
                                        output=outputfn, oformat='SQLite',
                                        co=['OVERWRITE=YES'],
                                        bandname=bandname)
                assert v.properties.getFeatureCount() == 11, \
                    'Error in geometry.extract() feature count (1)'

                v.io.close()
            else:
                v1 = pj.JimVect(outputfn)
                v2 = pj.geometry.extract(sampleid,
                                         jim, rule='mean',
                                         output='/vsimem/v2.sqlite',
                                         oformat='SQLite',
                                         co=['OVERWRITE=YES'],
                                         bandname=bandname)
                v = pj.geometry.join(v1, v2, output=outputfn, oformat='SQLite',
                                        co=['OVERWRITE=YES'], key=['fid'])
                assert v.properties.getFeatureCount() == 11, \
                    'Error in geometry.extract() feature count (2)'
                if band == 11:
                    assert len(v.properties.getFieldNames()) == nband + 2, \
                        'Error in geometry.extract() number of field names'
                v1.io.close()
                v2.io.close()
                v.io.close()
            jim.io.close()
        sampleid.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_extract_multiband():
        """Test the extractOgr method for jim with multibands."""
        jim0 = pj.Jim(rasterfn, band=[0, 1])
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()

        nband = jim0.properties.nrOfBand()

        planename = []
        for plane in range(0, jim0.properties.nrOfPlane()):
            planename.append('T' + str(plane))
        bandname = []
        for band in range(0, nband):
            bandname.append('B' + str(band))
        v = pj.geometry.extract(sampleid, jim0, rule='mean',
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'],
                                bandname=bandname)

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 11, \
            'Error in geometry.extract() feature count (1)'
        assert 'fid' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'fid' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert len(field_names) == nband + 2, \
            'Error in geometry.extract() field names (1)'

        newsamplefn = pj._get_random_path()
        sampleid.geometry.extract(jim0, rule='mean',
                                  output=newsamplefn, oformat='SQLite',
                                  co=['OVERWRITE=YES'],
                                  bandname=bandname)
        assert sampleid.properties.isEqual(v), \
            'Error in geometry.extract() function not equal to method'

        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)
        os.remove(newsamplefn)

    @staticmethod
    def test_extract_multiband_name():
        """Test the extractOgr method for jim with multibands with names."""
        jim0 = pj.Jim(rasterfn, band=[0, 1])
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()

        nband = jim0.properties.nrOfBand()

        planename = []
        for plane in range(0, jim0.properties.nrOfPlane()):
            planename.append('time' + str(plane))
        bandname = []
        jim0.properties.setDimension({'plane': planename, 'band': bandname})
        for band in range(0, nband):
            bandname.append('band' + str(band))
        v = pj.geometry.extract(sampleid, jim0, rule='mean',
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'])

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 11, \
            'Error in geometry.extract() feature count (1)'
        assert 'time0band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'time0band1' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'fid' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert len(field_names) == nband + 2, \
            'Error in geometry.extract() field names (1)'

        newsamplefn = pj._get_random_path()
        sampleid.geometry.extract(jim0, rule='mean',
                                  output=newsamplefn, oformat='SQLite',
                                  co=['OVERWRITE=YES'],
                                  bandname=bandname)
        assert sampleid.properties.isEqual(v), \
            'Error in geometry.extract() function not equal to method'

        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)
        os.remove(newsamplefn)

    @staticmethod
    def test_extract_plane():
        """Test the extractOgr method for multiplanes."""
        jim0 = pj.Jim(rasterfn, band= [0, 1], band2plane = True)
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid',
                                newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()

        nplane = jim0.properties.nrOfPlane()

        planename = []
        for plane in range(0, jim0.properties.nrOfPlane()):
            planename.append('T' + str(plane))
        bandname = []
        for band in range(0, jim0.properties.nrOfBand()):
            bandname.append('B' + str(band))

        jim0.properties.setDimension({'plane': planename, 'band': bandname})
        v = pj.geometry.extract(sampleid, jim0, rule='mean',
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'],
                                planename=planename,
                                bandname=bandname)

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 11, \
            'Error in geometry.extract() feature count (1)'
        assert 'fid' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'label' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert len(field_names) == nplane + 2, \
            'Error in geometry.extract() field names (1)'
        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_extract_plane_name():
        """Test the extractOgr method for multiplanes names."""
        jim0 = pj.Jim(rasterfn, band= [0, 1], band2plane = True)
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid',
                                newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()

        nplane = jim0.properties.nrOfPlane()

        planename = []
        for plane in range(0, jim0.properties.nrOfPlane()):
            planename.append('time' + str(plane))
        bandname = []
        for band in range(0, jim0.properties.nrOfBand()):
            bandname.append('band' + str(band))

        jim0.properties.setDimension({'plane': planename, 'band': bandname})
        v = pj.geometry.extract(sampleid, jim0, rule='mean',
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'])

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 11, \
            'Error in geometry.extract() feature count (1)'
        assert 'fid' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'label' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'time0band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'time1band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert len(field_names) == nplane + 2, \
            'Error in geometry.extract() field names (1)'
        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_extract_band_plane():
        """Test the extract method for multiband multiplanes."""
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()
        jim0 = pj.Jim(rasterfn, band=[0, 1])
        jim1 = pj.Jim(rasterfn, band=[6, 7])
        jim0.geometry.stackPlane(jim1)

        nband = jim0.properties.nrOfBand()
        nplane = jim0.properties.nrOfPlane()

        rule = ['mean', 'stdev']

        planename = []
        bandname = []
        for plane in range(0, nplane):
            planename.append('time' + str(plane))
        for band in range(0, nband):
            bandname.append('band' + str(band))
        v = pj.geometry.extract(sampleid, jim0, rule=['mean', 'stdev'],
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'],
                                bandname=bandname,
                                planename=planename)

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 11, \
            'Error in geometry.extract() feature count (1)'
        assert 'fid' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'label' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert len(field_names) == (nband + nplane) * len(rule) + 2, \
            'Error in geometry.extract() field names (1)'
        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_extract_band_plane_name():
        """Test the extract method for multiband multiplanes names."""
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()
        jim0 = pj.Jim(rasterfn, band=[0, 1])
        jim1 = pj.Jim(rasterfn, band=[6, 7])
        jim0.geometry.stackPlane(jim1)

        nband = jim0.properties.nrOfBand()
        nplane = jim0.properties.nrOfPlane()

        rule = ['mean', 'stdev']

        planename = []
        bandname = []
        for plane in range(0, nplane):
            planename.append('time' + str(plane))
        for band in range(0, nband):
            bandname.append('band' + str(band))
        jim0.properties.setDimension({'plane': planename, 'band': bandname})
        v = pj.geometry.extract(sampleid, jim0, rule=['mean', 'stdev'],
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'])

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 11, \
            'Error in geometry.extract() feature count (1)'
        assert 'fid' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'label' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'meantime0band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'meantime1band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'meantime0band1' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'meantime1band1' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'stdevtime0band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'stdevtime1band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'stdevtime0band1' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'stdevtime1band1' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert len(field_names) == (nband + nplane) * len(rule) + 2, \
            'Error in geometry.extract() field names (1)'
        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_extract_allpoints():
        """Test the extract method for allpoints."""
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid1',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()
        jim0 = pj.Jim(rasterfn, band=[0, 1])
        jim1 = pj.Jim(rasterfn, band=[6, 7])
        jim0.geometry.stackPlane(jim1)

        nband = jim0.properties.nrOfBand()
        nplane = jim0.properties.nrOfPlane()

        planename = []
        bandname = []
        for plane in range(0, jim0.properties.nrOfPlane()):
            planename.append('Time' + str(plane))
        for band in range(0, jim0.properties.nrOfBand()):
            bandname.append('Band' + str(band))
        v = pj.geometry.extract(sampleid, jim0, rule=['allpoints'],
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'],
                                copy='label',
                                bandname=bandname,
                                planename=planename)

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 307, \
            'Error in geometry.extract() feature count'
        assert len(field_names) == nband + nplane + 1, \
            'Error in geometry.extract() field names'
        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)

        sample = pj.JimVect(nutsfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid2',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()
        planename = []
        bandname = []
        for plane in range(0, jim0.properties.nrOfPlane()):
            planename.append('time' + str(plane))
        for band in range(0, jim0.properties.nrOfBand()):
            bandname.append('band' + str(band))
        v = pj.geometry.extract(sampleid, jim0, rule=['allpoints'],
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'],
                                copy='nuts3id',
                                classes=[255],
                                bandname=bandname,
                                planename=planename,
                                threshold=[10])

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 10, \
            'Error in geometry.extract() feature count'
        assert len(field_names) == nband + nplane + 1, \
            'Error in geometry.extract() field names'
        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_extract_allpoints_name():
        """Test the extract method for allpoints name."""
        sample = pj.JimVect(vectorfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid3',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()
        jim0 = pj.Jim(rasterfn, band=[0, 1])
        jim1 = pj.Jim(rasterfn, band=[6, 7])
        jim0.geometry.stackPlane(jim1)

        nband = jim0.properties.nrOfBand()
        nplane = jim0.properties.nrOfPlane()

        planename = []
        bandname = []
        for plane in range(0, jim0.properties.nrOfPlane()):
            planename.append('time' + str(plane))
        for band in range(0, jim0.properties.nrOfBand()):
            bandname.append('band' + str(band))
        jim0.properties.setDimension({'plane': planename, 'band': bandname})
        v = pj.geometry.extract(sampleid, jim0, rule=['allpoints'],
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'],
                                copy='label')

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 307, \
            'Error in geometry.extract() feature count'
        assert len(field_names) == nband + nplane + 1, \
            'Error in geometry.extract() field names'
        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)

        sample = pj.JimVect(nutsfn)
        sampleid = pj.JimVect(sample, output='/vsimem/sampleid4',
                              newfield='fid', co=['OVERWRITE=YES'])
        sample.io.close()
        v = pj.geometry.extract(sampleid, jim0, rule=['allpoints'],
                                output=outputfn, oformat='SQLite',
                                co=['OVERWRITE=YES'],
                                copy='nuts3id',
                                classes=[255],
                                threshold=[10])

        field_names = v.properties.getFieldNames()

        assert v.properties.getFeatureCount() == 10, \
            'Error in geometry.extract() feature count'
        assert 'time0band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'time1band0' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'time0band1' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert 'time1band1' in field_names, \
            'Error in geometry.extract() field names (1)'
        assert len(field_names) == nband + nplane + 1, \
            'Error in geometry.extract() field names'
        sampleid.io.close()
        v.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_extract_image():
        """Test the extract method for thematic reference."""
        jim0 = pj.Jim(rasterfn, band=[0, 1, 2, 3])
        jim0.geometry.warp('epsg:32632')
        reference = pj.Jim(clc, dx=500, dy=500)
        bbox = reference.properties.getBBox()
        jim0.geometry.crop(ulx=bbox[0], uly=bbox[1], lrx=bbox[2], lry=bbox[3])

        classes = [2, 12, 25, 41, 50]
        thresholds = [20, 20, 20, 20, 20]

        v = pj.geometry.extract(reference, jim0, srcnodata=[0],
                                rule='mean',
                                output=outputfn,
                                classes=classes,
                                threshold=thresholds,
                                bandname=['b02', 'b03', 'b04', 'b08'],
                                co=['OVERWRITE=YES'])

        assert v.properties.getFeatureCount() == 80, \
            'Error in extract method for thematic reference: feature count'
        assert 'label' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no label'
        assert 't0b02' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no t0b02'
        assert 't0b03' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no t0b03'
        assert 't0b04' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no t0b04'
        assert 't0b08' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no t0b08'
        v.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_extract_image_name():
        """Test the extract method for thematic reference name."""
        jim0 = pj.Jim(rasterfn, band=[0, 1, 2, 3])
        jim0.geometry.warp('epsg:32632')
        reference = pj.Jim(clc, dx=500, dy=500)
        bbox = reference.properties.getBBox()
        jim0.geometry.crop(ulx=bbox[0], uly=bbox[1], lrx=bbox[2], lry=bbox[3])

        classes = [2, 12, 25, 41, 50]
        thresholds = [20, 20, 20, 20, 20]

        bandname=['b02', 'b03', 'b04', 'b08']
        jim0.properties.setDimension({'plane': ['t0'], 'band': bandname})
        v = pj.geometry.extract(reference, jim0, srcnodata=[0],
                                rule='mean',
                                output=outputfn,
                                classes=classes,
                                threshold=thresholds,
                                co=['OVERWRITE=YES'])

        assert v.properties.getFeatureCount() == 80, \
            'Error in extract method for thematic reference: feature count'
        assert 'label' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no label'
        assert 't0b02' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no b2'
        assert 't0b03' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no b3'
        assert 't0b04' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no b4'
        assert 't0b08' in v.properties.getFieldNames(), \
            'Error in extract method for thematic reference: no b8'
        v.io.close()
        os.remove(outputfn)

    @staticmethod
    def test_crop():
        """Test crop...() functions and methods."""
        raster = pj.Jim(rasterfn, band=[0, 1])
        raster_stacked = pj.geometry.stackPlane(raster, raster)
        vector = pj.JimVect(vectorfn)
        raster_bbox = raster_stacked.properties.getBBox()
        raster_dx = raster_stacked.properties.getDeltaX()
        raster_dy = raster_stacked.properties.getDeltaY()
        vector_bbox = vector.properties.getBBox()

        # Test cropOgr()
        cropped = pj.geometry.cropOgr(raster_stacked, vector)
        raster_stacked.geometry.cropOgr(vector)
        raster_bbox_cropped = raster_stacked.properties.getBBox()

        mod_x = (raster_bbox_cropped[2] - raster_bbox_cropped[0]) % raster_dx
        mod_y = (raster_bbox_cropped[3] - raster_bbox_cropped[1]) % raster_dy

        assert raster_stacked.properties.isEqual(cropped), \
            'Inconsistency in geometry.cropOgr() ' \
            '(method returns different result than function)'
        assert raster_bbox != raster_bbox_cropped, \
            'Error in geometry.cropOgr() ' \
            '(BBox not changed after crop)'
        assert raster_bbox_cropped[:2] == vector_bbox[:2] and mod_x == 0 and \
               mod_y == 0, \
            'Error in geometry.cropOgr() ' \
            '(new BBox values not equal to the vector one)'

        # Test cropPlane()
        raster2 = pj.Jim(tiles[0])
        raster_stacked2 = pj.geometry.stackPlane(raster2, raster2, raster2)
        planes12 = pj.geometry.cropPlane(raster_stacked2, [0, 1])
        raster_stacked2.geometry.cropPlane([0, 1])

        assert planes12.properties.isEqual(raster_stacked2), \
            'Error in geometry.cropPlane() ' \
            '(method returns different result than function for list argument)'
        assert planes12.properties.nrOfPlane() == 2, \
            'Error in geometry.cropPlane() ' \
            '(not cropped to the right amount of planes for list argument)'

        plane1 = pj.geometry.cropPlane(raster_stacked2, 0)
        raster_stacked2.geometry.cropPlane(0)

        assert plane1.properties.isEqual(raster_stacked2), \
            'Error in geometry.cropPlane() ' \
            '(method returns different result than function for int argument)'
        assert plane1.properties.nrOfPlane() == 1, \
            'Error in geometry.cropPlane() ' \
            '(not cropped to the right amount of planes for int argument)'
        assert plane1.properties.isEqual(raster2), \
            'Error in geometry.cropPlane() ' \
            '(returned planes were changed)'

        #Test crop()
        raster = pj.Jim(rasterfn, band=[0, 1])
        raster_stacked = pj.geometry.stackPlane(raster, raster)
        bbox_orig = raster_stacked.properties.getBBox()
        bbox_1 = [bbox_orig[0] + raster_dx, bbox_orig[1] - raster_dy,
                bbox_orig[2] - raster_dx, bbox_orig[3] + raster_dy]

        raster_1 = pj.geometry.crop(raster_stacked, bbox = bbox_1)
        assert raster_1.properties.getBBox() == bbox_1, \
            'Error in geometry.crop() ' \
            'raster_1.properties.getBox() != bbox_1'
        assert raster_1.properties.nrOfCol() == raster.properties.nrOfCol() - 2, \
            'Error in geometry.crop() ' \
            'raster_1.properties.nrOfCol() != raster.properties.nrOfCol() - 2'
        assert raster_1.properties.nrOfRow() == raster.properties.nrOfRow() - 2, \
            'Error in geometry.crop() ' \
            'raster_1.properties.nrOfRow() != raster.properties.nrOfRow() - 2'
        assert raster_1.np(0)[0, 0, 0] == raster_stacked.np(0)[0, 1, 1], \
            'Error in geometry.crop() ' \
            'raster_1.np(0)[0, 0] == raster.np()[1, 1]'
        assert raster_1.np(0)[0, -1, -1] == raster_stacked.np(0)[0, -2, -2], \
            'Error in geometry.crop() ' \
            'raster_1.np(0)[-1, -1] == raster.np()[-2, -2]'
        assert raster_1.np(0)[1, 0, 0] == raster_stacked.np(0)[1, 1, 1], \
            'Error in geometry.crop() ' \
            'raster_1.np(0)[0, 0] == raster.np()[1, 1]'
        assert raster_1.np(0)[1, -1, -1] == raster_stacked.np(0)[1, -2, -2], \
            'Error in geometry.crop() ' \
            'raster_1.np(0)[-1, -1] == raster.np()[-2, -2]'
        assert raster_1.np(1)[0, 0, 0] == raster_stacked.np(1)[0, 1, 1], \
            'Error in geometry.crop() ' \
            'raster_1.np(1)[0, 0] == raster.np()[1, 1]'
        assert raster_1.np(1)[0, -1, -1] == raster_stacked.np(1)[0, -2, -2], \
            'Error in geometry.crop() ' \
            'raster_1.np(1)[-1, -1] == raster.np()[-2, -2]'
        assert raster_1.np(1)[1, 0, 0] == raster_stacked.np(1)[1, 1, 1], \
            'Error in geometry.crop() ' \
            'raster_1.np(1)[0, 0] == raster.np()[1, 1]'
        assert raster_1.np(1)[1, -1, -1] == raster_stacked.np(1)[1, -2, -2], \
            'Error in geometry.crop() ' \
            'raster_1.np(1)[-1, -1] == raster.np()[-2, -2]'

        bbox_05 = [bbox_orig[0] + raster_dx/2, bbox_orig[1] - raster_dy/2, bbox_orig[2] - raster_dx/2, bbox_orig[3] + raster_dy/2]
        raster_1 = pj.geometry.crop(raster_stacked, bbox = bbox_05)
        assert raster_1.properties.getBBox() == bbox_05, \
            'Error in geometry.crop() ' \
            'raster_1.properties.getBox() != bbox_05'
        assert raster_1.properties.nrOfCol() == raster.properties.nrOfCol() - 1, \
            'Error in geometry.crop() ' \
            'raster_1.properties.nrOfCol() != raster.properties.nrOfCol() - 1'
        assert raster_1.properties.nrOfRow() == raster.properties.nrOfRow() - 1, \
            'Error in geometry.crop() ' \
            'raster_1.properties.nrOfRow() != raster.properties.nrOfRow() - 1'

        bbox_15 = [bbox_orig[0] + raster_dx/2, bbox_orig[1] - raster_dy/2, bbox_orig[2], bbox_orig[3]]
        raster_1 = pj.geometry.crop(raster_stacked, bbox = bbox_15)
        assert raster_1.properties.getBBox()[0] == bbox_05[0], \
            'Error in geometry.crop() ' \
            'raster_1.properties.getBox()[0] == bbox_05[0]'
        assert raster_1.properties.getBBox()[1] == bbox_05[1], \
            'Error in geometry.crop() ' \
            'raster_1.properties.getBox()[1] == bbox_05[1]'
        assert raster_1.properties.getBBox()[2] == bbox_orig[2] + raster_dx/2, \
            'Error in geometry.crop() ' \
            'raster_1.properties.getBox()[2] == bbox_orig[2] + raster_dx/2'
        assert raster_1.properties.getBBox()[3] == bbox_orig[3] - raster_dy/2, \
            'Error in geometry.crop() ' \
            'raster_1.properties.getBox()[3] == bbox_orig[3] - raster_dy/2'
        assert raster_1.properties.nrOfCol() == raster.properties.nrOfCol(), \
            'Error in geometry.crop() ' \
            'raster_1.properties.nrOfCol() != raster.properties.nrOfCol()'
        assert raster_1.properties.nrOfRow() == raster.properties.nrOfRow(), \
            'Error in geometry.crop() ' \
            'raster_1.properties.nrOfRow() != raster.properties.nrOfRow()'

        #align
        raster_1 = pj.geometry.crop(raster_stacked, bbox = bbox_1, align = True)
        #should be alrearaster_dy aligned
        assert raster_1.properties.getBBox()[0] == bbox_1[0], \
            'Error in geometry.raster_1.properties.getBBox()() ' \
            'raster_1.properties.getBBox()[0] != bbox_1[0]'
        assert raster_1.properties.getBBox()[1] == bbox_1[1], \
            'Error in geometry.raster_1.properties.getBBox()() ' \
            'raster_1.properties.getBBox()[1] != bbox_1[1]'
        assert raster_1.properties.getBBox()[2] == bbox_1[2], \
            'Error in geometry.raster_1.properties.getBBox()() ' \
            'raster_1.properties.getBBox()[2] != bbox_1[2]'
        assert raster_1.properties.getBBox()[3] == bbox_1[3], \
            'Error in geometry.crop() ' \
            'raster_1.properties.getBBox()[3] != bbox_1[3]'

        raster_1 = pj.geometry.crop(raster_stacked, bbox = bbox_05, align = True)
        #should be original bbox
        assert raster_1.properties.getBBox()[0] == bbox_orig[0], \
            'Error in geometry.raster_1.properties.getBBox()() ' \
            'raster_1.properties.getBBox()[0] != bbox_orig[0]'
        assert raster_1.properties.getBBox()[1] == bbox_orig[1], \
            'Error in geometry.raster_1.properties.getBBox()() ' \
            'raster_1.properties.getBBox()[1] != bbox_orig[1]'
        assert raster_1.properties.getBBox()[2] == bbox_orig[2], \
            'Error in geometry.raster_1.properties.getBBox()() ' \
            'raster_1.properties.getBBox()[2] != bbox_orig[2]'
        assert raster_1.properties.getBBox()[3] == bbox_orig[3], \
            'Error in geometry.crop() ' \
            'raster_1.properties.getBBox()[3] != bbox_orig[3]'


        raster_15 = pj.geometry.crop(raster_stacked, bbox = bbox_15, align = True)
        #should be original bbox
        assert raster_15.properties.getBBox()[0] == bbox_orig[0], \
            'Error in geometry.raster_15.properties.getBBox()() ' \
            'raster_15.properties.getBBox()[0] != bbox_orig[0]'
        assert raster_15.properties.getBBox()[1] == bbox_orig[1], \
            'Error in geometry.raster_15.properties.getBBox()() ' \
            'raster_15.properties.getBBox()[1] != bbox_orig[1]'
        assert raster_15.properties.getBBox()[2] == bbox_orig[2], \
            'Error in geometry.raster_15.properties.getBBox()() ' \
            'raster_15.properties.getBBox()[2] != bbox_orig[2]'
        assert raster_15.properties.getBBox()[3] == bbox_orig[3], \
            'Error in geometry.crop() ' \
            'raster_15.properties.getBBox()[3] != bbox_orig[3]'
        raster_crop = pj.geometry.crop(raster, bbox = raster_bbox_cropped)
        assert raster_crop.properties.getBBox() == raster_bbox_cropped, \
            'Error in geometry.crop() ' \
            'raster_crop.properties.getBox() != raster_bbox_cropeed'
        raster_stacked2.properties.getBBox()  == raster_crop.properties.getBBox(), \
            'Error in geometry.crop() ' \
            'function crop != method crop'


    @staticmethod
    def test_coords_transformations():
        """Test geo2image() and image2geo() functions and methods."""
        jim = pj.Jim(tiles[0])
        geo_ulx, geo_uly, geo_lrx, geo_lry = jim.properties.getBBox()
        delta_x = jim.properties.getDeltaX()
        delta_y = jim.properties.getDeltaY()
        ulx_pix_center = geo_ulx + delta_x / 2
        uly_pix_center = geo_uly - delta_y / 2

        im_x, im_y = pj.geometry.geo2image(jim, geo_ulx, geo_uly)
        im_x_m, im_y_m = jim.geometry.geo2image(geo_ulx, geo_uly)

        assert im_x == im_x_m and im_y == im_y_m, \
            'Inconsistency in geometry.geo2image() ' \
            '(method returns different result than function)'
        assert im_x == 0 and im_y == 0, \
            'Error in geometry.geo2image()' \
            '(geo2image(BBox[0], BBox[1]) did not return zeros)'

        im_x, im_y = pj.geometry.geo2image(jim, 0, 0)
        im_x_m, im_y_m = jim.geometry.geo2image(0, 0)

        assert im_x == im_x_m and im_y == im_y_m, \
            'Inconsistency in geometry.geo2image() ' \
            '(method returns different result than function)'
        assert im_x == -geo_ulx / delta_x and im_y == geo_uly / delta_y, \
            'Error in geometry.geo2image()' \
            '(geo2image(0, 0) did not return original values divided by dX/dY)'

        geo_ulx_2, geo_uly_2 = pj.geometry.image2geo(jim, 0, 0)
        geo_ulx_2_m, geo_uly_2_m = jim.geometry.image2geo(0, 0)

        assert geo_ulx_2 == geo_ulx_2_m and geo_uly_2 == geo_uly_2_m, \
            'Inconsistency in geometry.geo2image() ' \
            '(method returns different result than function)'
        assert geo_ulx_2 == ulx_pix_center and geo_uly_2 == uly_pix_center, \
            'Error in geometry.image2geo()' \
            '(geo2image(0, 0) did not return original values divided by dX/dY)'

        geo_ulx_2, geo_uly_2 = pj.geometry.image2geo(jim, im_x, im_y)
        geo_ulx_2_m, geo_uly_2_m = jim.geometry.image2geo(im_x, im_y)

        assert geo_ulx_2 == geo_ulx_2_m and geo_uly_2 == geo_uly_2_m, \
            'Inconsistency in geometry.geo2image() ' \
            '(method returns different result than function)'
        assert geo_ulx_2 == 5 and geo_uly_2 == -5, \
            'Error in geometry.image2geo()' \
            '(geo2image(0, 0) did not return original values divided by dX/dY)'

    @staticmethod
    def test_image_frames():
        """Test imageFrame...() functions and methods."""
        nrow = ncol = 50
        nband = nplane = 2
        jim = pj.Jim(nrow=nrow, ncol=ncol, nband=nband, nplane=nplane,
                     otype='Byte', uniform=[0, 2], seed=0)

        # Test imageFrameAdd()
        #      (for 1-band Jim, see test below)

        added = pj.geometry.imageFrameAdd(jim, 1, 2, 1, 2, 1, 2, 5)
        jim.geometry.imageFrameAdd(1, 2, 1, 2, 1, 2, 5)

        assert jim.properties.isEqual(added), \
            'Inconsistency in geometry.imageFrameAdd() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of cols not raised or not raised to the right number)'
        assert jim.properties.nrOfRow() == nrow + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of rows not raised or not raised to the right number)'
        assert jim.properties.nrOfPlane() == nplane + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of planes not raised or not raised to the right number)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of bands changed)'
        assert jim[0, 25, 25].np() == 5, \
            'Error in geometry.imageFrameAdd() ' \
            '(value not used for new planes)'
        assert jim[1, 0, 0].np() == 5, \
            'Error in geometry.imageFrameAdd() ' \
            '(value not used for the frame)'
        assert 0 <= jim[1, 1, 1].np()[0] <= 1, \
            'Error in geometry.imageFrameAdd() ' \
            '(values in the original image changed)'

        # Test imageFrameSubtract()
        #      (for 1-band Jim, see test below)

        subtracted = pj.geometry.imageFrameSubtract(jim, 1, 2, 1, 2, 1, 2)
        jim.geometry.imageFrameSubtract(1, 2, 1, 2, 1, 2)

        assert jim.properties.isEqual(subtracted), \
            'Inconsistency in geometry.imageFrameSubtract() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of cols not raised or not raised to the right number)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of rows not raised or not raised to the right number)'
        assert jim.properties.nrOfPlane() == nplane, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of planes not raised or not raised to the right number)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of bands changed)'
        assert (jim.np() == added.np()[1:3, 1:-2, 1:-2]).all(), \
            'Error in geometry.imageFrameSubtract() ' \
            '(changed values in the original Jim)'

        # Test imageFrameSet()

        setted = pj.geometry.imageFrameSet(jim, 1, 2, 1, 2, 1, 0, 5)
        jim.geometry.imageFrameSet(1, 2, 1, 2, 1, 0, 5)

        assert jim.properties.isEqual(setted), \
            'Inconsistency in geometry.imageFrameSet() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.imageFrameSet() ' \
            '(number of cols changed)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.imageFrameSet() ' \
            '(number of rows changed)'
        assert jim.properties.nrOfPlane() == nplane, \
            'Error in geometry.imageFrameSet() ' \
            '(number of planes changed)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameSet() ' \
            '(number of bands changed)'
        assert jim[1, 0, 0].np() == 5, \
            'Error in geometry.imageFrameSet() ' \
            '(value not used for the frame)'
        assert 0 <= jim[1, 1, 1].np()[0] <= 1,\
            'Error in geometry.imageFrameSet() ' \
            '(values outside the frame changed)'

        # Test imageFrameSet() with a specified band

        setted = pj.geometry.imageFrameSet(jim, 1, 2, 1, 2, 0, 1, 10, band=0)
        jim.geometry.imageFrameSet(1, 2, 1, 2, 0, 1, 10, band=0)

        assert jim.properties.isEqual(setted), \
            'Inconsistency in geometry.imageFrameSet() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.imageFrameSet() ' \
            '(number of cols changed)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.imageFrameSet() ' \
            '(number of rows changed)'
        assert jim.properties.nrOfPlane() == nplane, \
            'Error in geometry.imageFrameSet() ' \
            '(number of planes changed)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameSet() ' \
            '(number of bands changed)'
        assert jim[1, 0, 0].np() == 10, \
            'Error in geometry.imageFrameSet() ' \
            '(value not used for the frame)'
        assert jim[0, 1, 1].np() == 5,\
            'Error in geometry.imageFrameSet() ' \
            '(values outside the frame changed)'

        # Test imageFrameAdd() for 1-band Jim

        jim.geometry.cropBand(0)
        nband = 1

        added = pj.geometry.imageFrameAdd(jim, 1, 2, 1, 2, 1, 2, 10)
        jim.geometry.imageFrameAdd(1, 2, 1, 2, 1, 2, 10)

        assert jim.properties.isEqual(added), \
            'Inconsistency in geometry.imageFrameAdd() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of cols not raised or not raised to the right number)'
        assert jim.properties.nrOfRow() == nrow + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of rows not raised or not raised to the right number)'
        assert jim.properties.nrOfPlane() == nplane + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of planes not raised or not raised to the right number)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of bands changed)'
        assert jim[0, 25, 25].np() == 10, \
            'Error in geometry.imageFrameAdd() ' \
            '(value not used for new planes)'
        assert jim[1, 0, 0].np() == 10, \
            'Error in geometry.imageFrameAdd() ' \
            '(value not used for the frame)'
        assert jim[1, 2, 2].np() == 5, \
            'Error in geometry.imageFrameAdd() ' \
            '(values in the original image changed)'

        # Test imageFrameSubtract() for 1-band Jim

        subtracted = pj.geometry.imageFrameSubtract(jim, 1, 2, 1, 2, 1, 2)
        jim.geometry.imageFrameSubtract(1, 2, 1, 2, 1, 2)

        assert jim.properties.isEqual(subtracted), \
            'Inconsistency in geometry.imageFrameSubtract() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of cols not raised or not raised to the right number)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of rows not raised or not raised to the right number)'
        assert jim.properties.nrOfPlane() == nplane, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of planes not raised or not raised to the right number)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of bands changed)'
        assert (jim.np() == added.np()[1:3, 1:-2, 1:-2]).all(), \
            'Error in geometry.imageFrameSubtract() ' \
            '(changed values in the original Jim)'

        jim = pj.Jim(ncol = 1024, nrow = 1024, nplane = 3, otype = 'GDT_Byte')
        jim.pixops.setData(1)
        jim1 = pj.geometry.imageFrameSet(jim, 1, 1, 1, 1, 0, 0, 0)

        jim[:,0:1,:] = 0
        jim[:,-1:,:] = 0
        jim[:,:,0:1] = 0
        jim[:,:,-1:] = 0

        assert(jim.properties.isEqual(jim1)),\
            'Error in geometry.imageFrameSet() ' \
            '(result different from __setitem__)'

    # def test_image_inserts():
    #     """Test imageFrame...() functions and methods."""
    #     nrow = ncol = 50
    #     nband = nplane = 2
    #     jim = pj.Jim(nrow=nrow, ncol=ncol, nband=nband, nplane=nplane,
    #                  otype='Byte', uniform=[0, 2])
    #     jim_copy = pj.Jim(jim)
    #
    #     ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
    #     ngb[0, 1] = 1
    #     ngb[1, 0] = 1
    #     ngb[1, 2] = 1
    #     ngb[2, 1] = 1
    #
    #     # Test imageInsert()
    #     #      (for 1-band Jim, see test below)
    #
    #     inserted = pj.geometry.imageInsert(jim, ngb, 1, 1, 0)
    #     jim.geometry.imageInsert(ngb, 1, 1, 0)

    @staticmethod
    def test_plotLine():
        """Test the plotLine() function and method."""
        nrow = ncol = 10
        nband = 2
        jim = pj.Jim(nrow=nrow, ncol=ncol, nband=nband,
                     otype='Byte', uniform=[0, 50], seed=0)

        avg = jim.stats.getStats('mean')['mean']

        plotted = pj.geometry.plotLine(jim, 0, 0, 9, 9, 100)
        jim.geometry.plotLine(0, 0, 9, 9, 100)

        new_avg = jim.stats.getStats('mean')['mean']

        assert jim.properties.isEqual(plotted), \
            'Inconsistency in geometry.plotLine() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.plotLine() ' \
            '(number of cols changed)'
        assert new_avg > avg, \
            'Error in geometry.plotLine() ' \
            '(average not higher than of the original image, but higher ' \
            'values were set for the line)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.plotLine() ' \
            '(number of rows changed)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.plotLine() ' \
            '(number of planes changed)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.plotLine() ' \
            '(number of bands changed)'
        assert (jim.np()[range(5), range(5)] == 100).all(), \
            'Error in geometry.plotLine() ' \
            '(values not changed to the right value)'
        assert (jim.np()[range(1, 5), range(4)] != 100).all(), \
            'Error in geometry.plotLine() ' \
            '(values outside the line changed)'

        # Test wrong calls
        jim = pj.Jim(nrow=nrow, ncol=ncol, nband=nband, nplane=2,
                     otype='Byte', uniform=[0, 50], seed=0)

        try:
            _ = pj.geometry.plotLine(jim, 0, 0, 9, 9, 100)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.plotLine() ' \
            'function where the Jim argument is a multi-plane object'

        # Test wrong calls
        try:
            jim.geometry.plotLine(0, 0, 9, 9, 100)
            raised = False
        except pj.exceptions.JimInnerParametersError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.plotLine() ' \
            'method where the Jim argument is a multi-plane object'

    @staticmethod
    def test_polygonize():
        """Test the polygonize() function and method."""
        jim = pj.Jim(tiles[0])

        sub = int(jim.properties.nrOfCol() / 2 - 3)

        jim.geometry.imageFrameSubtract(sub, sub, sub, sub)

        jim[0, 0] = 1
        jim[0, 1] = 1

        pol1fn = pj._get_random_path()
        pol2fn = pj._get_random_path()
        pol1 = pj.geometry.polygonize(jim, pol1fn)
        pol2 = jim.geometry.polygonize(pol2fn)

        feature_count_func = pol1.properties.getFeatureCount()
        feature_count_meth = pol2.properties.getFeatureCount()

        nr_of_cells = jim.properties.nrOfCol() * jim.properties.nrOfRow()

        assert feature_count_func == feature_count_meth, \
            'Inconsistency in geometry.polygonize() ' \
            '(method returns different result than function)'
        assert pol1.properties.getBBox() == pol2.properties.getBBox() == \
               jim.properties.getBBox(), \
            'Error in geometry.polygonize() ' \
            '(BBox changed)'
        assert feature_count_func < nr_of_cells, \
            'Error in geometry.polygonize() ' \
            '(not less features in polygons than cells in raster)'

        os.remove(pol1fn)
        os.remove(pol2fn)

        # Test with the mask parameter
        mask = pj.Jim(jim, copy_data=False)
        mask[0, 0] = 1
        mask[0, 1] = 1
        mask[2, 2] = 1

        pol1_mask = pj.geometry.polygonize(jim, pol1fn,
                                           mask=mask)
        pol2_mask = jim.geometry.polygonize(pol2fn, mask=mask)

        feature_count_func_mask = pol1_mask.properties.getFeatureCount()
        feature_count_meth_mask = pol2_mask.properties.getFeatureCount()

        assert feature_count_func_mask == feature_count_meth_mask, \
            'Inconsistency in geometry.polygonize() ' \
            '(method returns different result than function when mask ' \
            'argument used)'
        assert pol1_mask.properties.getBBox() == \
               pol2_mask.properties.getBBox() == \
               jim[:3, :3].properties.getBBox(), \
            'Error in geometry.polygonize() ' \
            '(BBox not changed or changed wrongly when mask argument used)'
        assert 2 <= feature_count_func_mask < feature_count_func, \
            'Error in geometry.polygonize() ' \
            '(not less features in polygons when mask argument used)'

        os.remove(pol1fn)
        os.remove(pol2fn)

        # Test wrong calls
        try:
            _ = pj.geometry.polygonize(1, pj._get_random_path())
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.polygonize(jim, path) ' \
            'function where the jim argument is not an instance of a Jim ' \
            'object'

        try:
            _ = pj.geometry.polygonize(jim, pj._get_random_path(), mask=5)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.polygonize(jim, path, ' \
            'mask) function where the mask argument is not an instance of a ' \
            'Jim object'

        try:
            _ = jim.geometry.polygonize(pj._get_random_path(), mask='spam')
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.polygonize(path, mask) ' \
            'method where the mask argument is not an instance of a Jim object'

    @staticmethod
    def test_rasterize():
        """Test the rasterize() function and method for bytes."""
        jim0 = pj.Jim(rasterfn, band=0)
        sample = pj.JimVect(vectorfn)
        mask = pj.Jim(jim0, copy_data=False)
        mask.pixops.convert('GDT_Byte')
        mask.geometry.rasterize(sample, eo=['ATTRIBUTE=label'])

        rasterized = pj.geometry.rasterize(jim0, sample,
                                           eo=['ATTRIBUTE=label'])

        assert rasterized.properties.isEqual(mask), \
            'Error in geometry.rasterize() ' \
            '(function is not equal to method)'

        minmax = mask.stats.getStats(function=['min', 'max'])

        assert minmax['min'] == 0,\
            'Error in geometry.rasterize() min != 0'
        assert minmax['max'] == 2,\
            'Error in geometry.rasterize() max != 2'

        # Test the rasterize() function and method for double.
        jim0 = pj.Jim(rasterfn, band=0)
        jim0.pixops.convert('GDT_Float64')
        sample = pj.JimVect(vectorfn)
        mask = pj.Jim(jim0, copy_data=False)
        mask.geometry.rasterize(sample, eo=['ATTRIBUTE=label'])

        rasterized = pj.geometry.rasterize(jim0, sample,
                                           eo=['ATTRIBUTE=label'])

        assert rasterized.properties.isEqual(mask), \
            'Error in geometry.rasterize() ' \
            '(function is not equal to method)'

        minmax = mask.stats.getStats(function=['min', 'max'])

        assert minmax['min'] == 0,\
            'Error in geometry.rasterize() min != 0'
        assert minmax['max'] == 2,\
            'Error in geometry.rasterize() max != 2'

        # Test with ln parameter
        rasterized2 = pj.geometry.rasterize(jim0, sample,
                                            eo=['ATTRIBUTE=label'],
                                            ln='training')

        assert rasterized.properties.isEqual(rasterized2), \
            'Error in geometry.rasterize() ' \
            '(when using "ln" parameter to get all layers, the object is ' \
            'different than with ln=None)'

        jim0.geometry.rasterize(sample, eo=['ATTRIBUTE=label'], ln='training')

        assert rasterized.properties.isEqual(jim0), \
            'Inconsistency in geometry.rasterize() ' \
            '(method returns different result than function)'

        # Test wrong calls

        try:
            _ = pj.geometry.rasterize('spam', sample)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.rasterize(jim, jimvect) ' \
            'function where the jim argument is not an instance of a Jim ' \
            'object'

        try:
            _ = pj.geometry.rasterize(jim0, 'spam')
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.rasterize(jim, jimvect) ' \
            'function where the jimvect argument is not an instance of a ' \
            'JimVect object'

        try:
            _ = jim0.geometry.rasterize('spam')
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.rasterize(jim, jimvect) ' \
            'function where the jimvect argument is not an instance of a ' \
            'JimVect object'

        sample.io.close()
        jim0.io.close()
        mask.io.close()

    @staticmethod
    def test_reducePlane():
        """Test the reducePlane() function and method."""
        nr_of_row = nr_of_col = 10
        min = 0
        max = 10
        nodata = int((max + min) / 2)

        # Test with no rule specified (=overwrite)
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3, otype='Byte',
                     uniform=[min, max])
        jim_copy = pj.Jim(jim)
        last_plane = pj.Jim(jim[-1])

        reduced = pj.geometry.reducePlane(jim)
        reduced2 = pj.geometry.reducePlane(jim, rule='overwrite')
        jim.geometry.reducePlane()
        jim_copy.geometry.reducePlane(rule='overwrite')

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function when no rule ' \
            'specified)'
        assert jim_copy.properties.isEqual(reduced2), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function for ' \
            'rule="overwrite")'
        assert jim.properties.isEqual(reduced2), \
            'Error in geometry.reducePlane() ' \
            '(different result for no rule and rule="overwrite" - the ' \
            'default value not passed)'
        assert jim.properties.isEqual(last_plane), \
            'Error in geometry.reducePlane() ' \
            '(returned Jim is not equal to the last plane for ' \
            'rule="overwrite")'

        # Test with sum rule via call back function
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nband = 3, nplane=1, otype='UInt16')
        jim.pixops.setData(1, bands=[0])
        jim.pixops.setData(2, bands=[1])
        jim.pixops.setData(3, bands=[2])
        jim.geometry.stackPlane(pj.Jim(jim))

        jim1 = pj.geometry.cropPlane(jim, 0)
        for band in range(0, jim.properties.nrOfBand()):
            jim1.np(band)[:] = np.sum(jim.np(band), axis = 0)

        def getSum(reduced, plane):
            reduced+=plane
            return reduced

        jim.geometry.reducePlane(getSum)

        assert jim.properties.isEqual(jim1) ,\
            'Inconsistency in geometry.reducePlane() with sum callback, ' \
            'not equal to Numpy implementation'
        assert jim.stats.getStats(['min','max']) == {'min': [2,4,6], 'max': [2,4,6]} ,\
            'Inconsistency in geometry.reducePlane() with sum callback'

        # Test with rule == 'max'
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3, otype='Byte',
                     uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, 'max')
        jim.geometry.reducePlane('max')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == max, \
            'Error in geometry.reducePlane() ' \
            '(rule="max" did not return max value for all the planes)'
        assert stats_reduced['min'] >= stats['min'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the minimum value of returned object is not >=' \
            ' the minimum of the original object)'
        assert stats_reduced['max'] == stats['max'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the maximum value of returned object is not ' \
            'equal to the maximum of the original object)'
        assert stats_reduced['mean'] > stats['mean'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the mean value of returned object is not > ' \
            'the mean of the original object)'

        # Test with rule == 'max' and multibands
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, 'max')
        jim.geometry.reducePlane('max')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == max, \
            'Error in geometry.reducePlane() ' \
            '(rule="max" did not return max value for all the planes)'
        assert all(stats_reduced['min'][i] >= stats['min'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the minimum value of returned object is not >=' \
            ' the minimum of the original object)'
        assert all(stats_reduced['max'][i] == stats['max'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the maximum value of returned object is not ' \
            'equal to the maximum of the original object)'
        assert all(stats_reduced['mean'][i] > stats['mean'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the mean value of returned object is not > ' \
            'the mean of the original object)'

        # Test with rule == 'max' and band specified
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        jim.np(1)[2, 0, 0] = int(2.5 * max)
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, 'max', ref_band=0)
        jim.geometry.reducePlane('max', ref_band=0)

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == max, \
            'Error in geometry.reducePlane() ' \
            '(rule="max" did not return max value for all the planes)'
        assert jim.np(1)[0, 0] == 2 * max, \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="max", the used indices of max values are not he ones' \
            ' from the ref_band)'
        assert stats_reduced['min'][0] >= stats['min'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the minimum value of returned object is not >=' \
            ' the minimum of the original object)'
        assert stats_reduced['mean'][0] > stats['mean'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the mean value of returned object is not > ' \
            'the mean of the original object)'

        # Test with rule == 'max' and nodata and band specified
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim[:, 2, 2] = nodata
        jim.np(1)[:, 2, 2] = nodata + 1
        jim[0, 3, 3] = nodata
        jim[1, 3, 3] = max
        jim[2, 3, 3] = nodata + 1
        jim.np(1)[0, 3, 3] = max
        jim.np(1)[1, 3, 3] = nodata
        jim.np(1)[2, 3, 3] = nodata + 1
        jim.np(0)[0, 4, 4] = nodata + 1
        jim.np(1)[0, 4, 4] = nodata
        jim.np(1)[1, 4, 4] = nodata
        jim.np(1)[2, 4, 4] = nodata
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, 'max', ref_band=0,
                                          nodata=nodata)
        jim.geometry.reducePlane('max', ref_band=0, nodata=nodata)

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == max, \
            'Error in geometry.reducePlane() ' \
            '(rule="max" did not return max value for all the planes)'
        assert jim.np(1)[0, 0] == 2 * max, \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="max", the used indices of max values are not he ones' \
            ' from the ref_band)'
        assert jim.np(0)[3, 3] == max, \
            'Error in geometry.reducePlane(nodata) ' \
            '(not ignoring the nodata values)'
        assert jim.np(1)[3, 3] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(not ignoring the values with indices of nodata from the ' \
            'ref_band for other bands)'
        assert jim.np(1)[4, 4] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(ignoring also the nodata values in other bands than ref_band)'
        assert jim.np(0)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata) ' \
            '(the returned object not containing nodata on a place where ' \
            'nodata was in all planes)'
        assert jim.np(1)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(the returned object not containing nodata in all bands on a ' \
            'place where nodata was in all planes in the ref_band)'
        assert stats_reduced['min'][0] >= stats['min'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the minimum value of returned object is not >=' \
            ' the minimum of the original object)'
        assert stats_reduced['mean'][0] > stats['mean'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the mean value of returned object is not > ' \
            'the mean of the original object)'

        # Test with rule == 'min'
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3, otype='Byte',
                     uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, rule='min')
        jim.geometry.reducePlane(rule='min')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == min, \
            'Error in geometry.reducePlane() ' \
            '(rule="min" did not return max value for all the planes)'
        assert stats_reduced['max'] <= stats['max'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the maximum value of returned object is not <=' \
            ' the maximum of the original object)'
        assert stats_reduced['min'] == stats['min'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the minimum value of returned object is not ' \
            'equal to the minimum of the original object)'
        assert stats_reduced['mean'] < stats['mean'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the mean value of returned object is not < ' \
            'the mean of the original object)'

        # Test with rule == 'min' and multibands
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, rule='min')
        jim.geometry.reducePlane(rule='min')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == min, \
            'Error in geometry.reducePlane() ' \
            '(rule="min" did not return min value for all the planes)'
        assert all(stats_reduced['max'][i] <= stats['max'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the maximum value of returned object is not <=' \
            ' the maximum of the original object)'
        assert all(stats_reduced['min'][i] == stats['min'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the minimum value of returned object is not ' \
            'equal to the minimum of the original object)'
        assert all(stats_reduced['mean'][i] < stats['mean'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="max", the mean value of returned object is not < ' \
            'the mean of the original object)'

        # Test with rule == 'min' and band specified
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = nodata
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, ref_band=0, rule='min')
        jim.geometry.reducePlane(ref_band=0, rule='min')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == min, \
            'Error in geometry.reducePlane() ' \
            '(rule="min" did not return min value for all the planes)'
        assert jim.np(1)[0, 0] == 3 * max, \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="min", the used indices of max values are not he ones' \
            ' from the ref_band)'
        assert stats_reduced['max'][0] <= stats['max'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the maximum value of returned object is not <=' \
            ' the maximum of the original object)'
        assert stats_reduced['mean'][0] < stats['mean'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the mean value of returned object is not < ' \
            'the mean of the original object)'

        # Test with rule == 'min' and nodata and ref_band specified
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim[:, 2, 2] = nodata
        jim.np(1)[:, 2, 2] = nodata + 1
        jim[0, 3, 3] = nodata
        jim[1, 3, 3] = min
        jim[2, 3, 3] = nodata - 1
        jim.np(1)[0, 3, 3] = min
        jim.np(1)[1, 3, 3] = nodata
        jim.np(1)[2, 3, 3] = nodata - 1
        jim.np(0)[0, 4, 4] = nodata + 1
        jim.np(1)[0, 4, 4] = nodata
        jim.np(1)[1, 4, 4] = nodata
        jim.np(1)[2, 4, 4] = nodata
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(
            jim, ref_band=0, nodata=nodata, rule='min')
        jim.geometry.reducePlane(ref_band=0, nodata=nodata, rule='min')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == min, \
            'Error in geometry.reducePlane() ' \
            '(rule="max" did not return max value for all the planes)'
        assert jim.np(1)[0, 0] == 3 * max, \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="min", the used indices of max values are not he ones' \
            ' from the ref_band)'
        assert jim.np(0)[3, 3] == min, \
            'Error in geometry.reducePlane(nodata) ' \
            '(not ignoring the nodata values)'
        assert jim.np(1)[3, 3] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(not ignoring the values with indices of nodata from the ' \
            'ref_band for other bands)'
        assert jim.np(1)[4, 4] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(ignoring also the nodata values in other bands than ref_band)'
        assert jim.np(0)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata) ' \
            '(the returned object not containing nodata on a place where ' \
            'nodata was in all planes)'
        assert jim.np(1)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(the returned object not containing nodata in all bands on a ' \
            'place where nodata was in all planes in the ref_band)'
        assert stats_reduced['max'][0] <= stats['max'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the maximum value of returned object is not <=' \
            ' the maximum of the original object)'
        assert stats_reduced['mean'][0] < stats['mean'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="min", the mean value of returned object is not < ' \
            'the mean of the original object)'

        # Test with rule == 'mean'
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3, otype='Byte',
                     uniform=[min, max])
        jim[2, 0, 0] = nodata
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, rule='mean')
        jim.geometry.reducePlane(rule='mean')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == (max + min + nodata) / 3, \
            'Error in geometry.reducePlane() ' \
            '(rule="mean" did not return mean value for all the planes)'
        assert stats_reduced['max'] <= stats['max'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the maximum value of returned object is not ' \
            '<= the maximum of the original object)'
        assert stats_reduced['min'] >= stats['min'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the minimum value of returned object is not ' \
            '>= the minimum of the original object)'

        # Test with rule == 'mean' and multibands
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = nodata
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, rule='mean')
        jim.geometry.reducePlane(rule='mean')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == int((max + min + nodata) / 3), \
            'Error in geometry.reducePlane() ' \
            '(rule="mean" did not return mean value for all the planes)'
        assert all(stats_reduced['max'][i] <= stats['max'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the maximum value of returned object is not ' \
            '<= the maximum of the original object)'
        assert all(stats_reduced['min'][i] >= stats['min'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the minimum value of returned object is not ' \
            '>= the minimum of the original object)'

        # Test with rule == 'mean' and band specified (2 planes)
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=2,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, ref_band=0, rule='mean')
        jim.geometry.reducePlane(ref_band=0, rule='mean')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == (max + min) / 2, \
            'Error in geometry.reducePlane() ' \
            '(rule="mean" did not return mean value for all the planes)'
        assert jim.np(1)[0, 0] == 5 * max / 2, \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="mean", the used indices of max values are not the ' \
            'ones from the ref_band)'
        assert all(stats_reduced['max'][i] <= stats['max'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the maximum value of returned object is ' \
            'not <= the maximum of the original object)'
        assert all(stats_reduced['min'][i] >= stats['min'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the min value of returned object is not >= ' \
            'the min of the original object)'

        # Test with rule == 'mean' and band specified (>2 planes)
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = (max + min) / 4
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        jim.np(1)[2, 0, 0] = 2 * max
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, ref_band=0, rule='mean')
        jim.geometry.reducePlane(ref_band=0, rule='mean')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function for more than ' \
            '2-planes Jim)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1) for more than 2-planes Jim'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed) for more than 2-planes Jim'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed for more than ' \
            '2-planes Jim)'
        assert jim.np()[0, 0] == int((max + min + int((max + min) / 4)) / 3), \
            'Error in geometry.reducePlane() ' \
            '(rule="mean" did not return mean value for all the planes for ' \
            'more than 2-planes Jim)'
        assert jim.np(1)[0, 0] == int(7 * max / 3), \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="mean", the used indices of max values are not the ' \
            'ones from the ref_band for more than 2-planes Jim)'
        assert all(stats_reduced['max'][i] <= stats['max'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the maximum value of returned object is not ' \
            '<= the maximum of the original object for more than 2-planes Jim)'
        assert all(stats_reduced['min'][i] >= stats['min'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the min value of returned object is not >= ' \
            'the min of the original object for more than 2-planes Jim)'

        # Test with rule == 'mean' and nodata and ref_band specified
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = nodata + 1
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim[:, 2, 2] = nodata
        jim.np(1)[:, 2, 2] = nodata + 1
        jim[0, 3, 3] = nodata
        jim[1, 3, 3] = nodata + 1
        jim[2, 3, 3] = max
        jim.np(1)[0, 3, 3] = 5 * max
        jim.np(1)[1, 3, 3] = nodata
        jim.np(1)[2, 3, 3] = max
        jim.np(0)[0, 4, 4] = nodata + 1
        jim.np(1)[0, 4, 4] = nodata
        jim.np(1)[1, 4, 4] = nodata
        jim.np(1)[2, 4, 4] = nodata
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        jim.np(1)[2, 0, 0] = nodata + 1
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(
            jim, ref_band=0, nodata=nodata, rule='mean')
        jim.geometry.reducePlane(ref_band=0, nodata=nodata, rule='mean')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == nodata, \
            'Error in geometry.reducePlane() ' \
            '(rule="mean" did not return mean value for all the planes)'
        assert jim.np(1)[0, 0] == int((5 * max + nodata + 1) / 3), \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="mean", the used indices of values are not ' \
            'the ones from the ref_band)'
        assert jim.np(0)[3, 3] == int(nodata + 1 + max) / 2, \
            'Error in geometry.reducePlane(nodata) ' \
            '(not ignoring the nodata values)'
        assert jim.np(1)[3, 3] == int((nodata + max) / 2), \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(not ignoring the values with indices of nodata from the ' \
            'ref_band for other bands)'
        assert jim.np(1)[4, 4] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(ignoring also the nodata values in other bands than ref_band)'
        assert jim.np(0)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata) ' \
            '(the returned object not containing nodata on a place where ' \
            'nodata was in all planes)'
        assert jim.np(1)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(the returned object not containing nodata in all bands on a ' \
            'place where nodata was in all planes in the ref_band)'
        assert stats_reduced['max'][0] <= stats['max'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the maximum value of returned object is ' \
            'not <= the maximum of the original object)'
        assert stats_reduced['min'][0] >= stats['min'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the mean value of returned object is not >= ' \
            'the mean of the original object)'

        # Test with rule == 'mean' and nodata specified
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = max
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim[:, 2, 2] = nodata
        jim.np(1)[:, 2, 2] = nodata + 1
        jim[0, 3, 3] = nodata
        jim[1, 3, 3] = nodata + 1
        jim[2, 3, 3] = nodata + 2
        jim.np(1)[0, 3, 3] = nodata + 1
        jim.np(1)[1, 3, 3] = nodata
        jim.np(1)[2, 3, 3] = nodata + 1
        jim.np(0)[0, 4, 4] = nodata + 1
        jim.np(1)[0, 4, 4] = nodata
        jim.np(1)[1, 4, 4] = nodata
        jim.np(1)[2, 4, 4] = nodata
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        jim.np(1)[2, 0, 0] = nodata
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, nodata=nodata, rule='mean')
        jim.geometry.reducePlane(nodata=nodata, rule='mean')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == int((2 * max + min) / 3), \
            'Error in geometry.reducePlane() ' \
            '(rule="mean" did not return mean value for all the planes)'
        assert jim.np(1)[0, 0] == 5 * max / 2, \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="mean", the used indices of max values are not ' \
            'the ones from the ref_band)'
        assert jim.np(0)[3, 3] == int((nodata + 1 + nodata + 2) / 2), \
            'Error in geometry.reducePlane(nodata) ' \
            '(not ignoring the nodata values)'
        assert jim.np(1)[3, 3] == nodata + 1, \
            'Error in geometry.reducePlane(nodata) ' \
            '(not ignoring the nodata values for other bands than band 0)'
        assert jim.np(1)[4, 4] == nodata, \
            'Error in geometry.reducePlane(nodata) ' \
            '(ignoring also the nodata values in other bands than ref_band)'
        assert jim.np(0)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata) ' \
            '(the returned object not containing nodata on a place where ' \
            'nodata was in all planes)'
        assert jim.np(1)[2, 2] == nodata + 1, \
            'Error in geometry.reducePlane(nodata) ' \
            '(the returned object not containing nodata in all bands on a ' \
            'place where nodata was in all planes in the ref_band)'
        assert stats_reduced['max'][0] <= stats['max'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the maximum value of returned object is not ' \
            '<= the maximum of the original object)'
        assert stats_reduced['mean'][0] < stats['mean'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="mean", the mean value of returned object is not < ' \
            'the mean of the original object)'

        # Test with rule == 'median'
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3, otype='Byte',
                     uniform=[min, max])
        jim[2, 0, 0] = nodata
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, rule='median')
        jim.geometry.reducePlane(rule='median')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == np.median([min, max, nodata]), \
            'Error in geometry.reducePlane() ' \
            '(rule="median" did not return median value for all the planes)'
        assert stats_reduced['max'] <= stats['max'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the maximum value of returned object is not' \
            ' <= the maximum of the original object)'
        assert stats_reduced['min'] >= stats['min'], \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the minimum value of returned object is not' \
            ' >= the minimum of the original object)'

        # Test with rule == 'median' and multibands
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = nodata
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, rule='median')
        jim.geometry.reducePlane(rule='median')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == np.median([min, max, nodata]), \
            'Error in geometry.reducePlane() ' \
            '(rule="median" did not return median value for all the planes)'
        assert all(stats_reduced['max'][i] <= stats['max'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the maximum value of returned object is not' \
            ' <= the maximum of the original object)'
        assert all(stats_reduced['min'][i] >= stats['min'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the minimum value of returned object is not' \
            ' >= the minimum of the original object)'

        # Test with rule == 'median' and band specified (2 planes)
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=2,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, ref_band=0, rule='median')
        jim.geometry.reducePlane(ref_band=0, rule='median')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == np.median([min, max]), \
            'Error in geometry.reducePlane() ' \
            '(rule="median" did not return median value for all the planes)'
        assert jim.np(1)[0, 0] == np.median([3 * max, 2 * max]), \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="median", the used indices of max values are not the ' \
            'ones from the ref_band)'
        assert all(stats_reduced['max'][i] <= stats['max'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the maximum value of returned object is ' \
            'not <= the maximum of the original object)'
        assert all(stats_reduced['min'][i] >= stats['min'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the min value of returned object is not >= ' \
            'the min of the original object)'

        # Test with rule == 'median' and band specified (>2 planes)
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = (max + min) / 4
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        jim.np(1)[2, 0, 0] = 2 * max
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, ref_band=0, rule='median')
        jim.geometry.reducePlane(ref_band=0, rule='median')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function for more than ' \
            '2-planes Jim)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1) for more than 2-planes Jim'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed) for more than 2-planes Jim'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed for more than ' \
            '2-planes Jim)'
        assert jim.np()[0, 0] == np.median([max, min, int((max + min) / 4)]), \
            'Error in geometry.reducePlane() ' \
            '(rule="median" did not return median value for all the planes' \
            ' for more than 2-planes Jim)'
        assert jim.np(1)[0, 0] == np.median([3 * max, 2 * max, 2 * max]), \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="median", the used indices of max values are not the ' \
            'ones from the ref_band for more than 2-planes Jim)'
        assert all(stats_reduced['max'][i] <= stats['max'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the maximum value of returned object is not' \
            ' <= the maximum of the original object for more than 2-planes ' \
            'Jim)'
        assert all(stats_reduced['min'][i] >= stats['min'][i] for i in range(
            len(stats['max']))), \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the min value of returned object is not >= ' \
            'the min of the original object for more than 2-planes Jim)'

        # Test with rule == 'median' and nodata and ref_band specified
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = nodata + 1
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim[:, 2, 2] = nodata
        jim.np(1)[:, 2, 2] = nodata + 1
        jim[0, 3, 3] = nodata
        jim[1, 3, 3] = nodata + 1
        jim[2, 3, 3] = max
        jim.np(1)[0, 3, 3] = 5 * max
        jim.np(1)[1, 3, 3] = nodata
        jim.np(1)[2, 3, 3] = max
        jim.np(0)[0, 4, 4] = nodata + 1
        jim.np(1)[0, 4, 4] = nodata
        jim.np(1)[1, 4, 4] = nodata
        jim.np(1)[2, 4, 4] = nodata
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        jim.np(1)[2, 0, 0] = nodata + 1
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(
            jim, ref_band=0, nodata=nodata, rule='median')
        jim.geometry.reducePlane(ref_band=0, nodata=nodata, rule='median')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == np.median([nodata + 1, min, max]), \
            'Error in geometry.reducePlane() ' \
            '(rule="median" did not return median value for all the planes)'
        assert jim.np(1)[0, 0] == \
               int(np.median([3 * max, 2 * max, nodata + 1])), \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="median", the used indices of values are not ' \
            'the ones from the ref_band)'
        assert jim.np(0)[3, 3] == int(np.median([nodata + 1, max])), \
            'Error in geometry.reducePlane(nodata) ' \
            '(not ignoring the nodata values)'
        assert jim.np(1)[3, 3] == int(np.median([nodata, max])), \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(not ignoring the values with indices of nodata from the ' \
            'ref_band for other bands for rule="median")'
        assert jim.np(1)[4, 4] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(ignoring also the nodata values in other bands than ref_band ' \
            'for rule="median")'
        assert jim.np(0)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata) ' \
            '(the returned object not containing nodata on a place where ' \
            'nodata was in all planes for rule="median")'
        assert jim.np(1)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata, ref_band) ' \
            '(the returned object not containing nodata in all bands on a ' \
            'place where nodata was in all planes in the ref_band for ' \
            'rule="median")'
        assert stats_reduced['max'][0] <= stats['max'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the maximum value of returned object is ' \
            'not <= the maximum of the original object)'
        assert stats_reduced['min'][0] >= stats['min'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the mean value of returned object is not >=' \
            ' the mean of the original object)'

        # Test with rule == 'median' and nodata specified
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='Byte', uniform=[min, max])
        jim[2, 0, 0] = max
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim[:, 2, 2] = nodata
        jim.np(1)[:, 2, 2] = nodata + 1
        jim[0, 3, 3] = nodata
        jim[1, 3, 3] = nodata + 1
        jim[2, 3, 3] = nodata + 2
        jim.np(1)[0, 3, 3] = nodata + 1
        jim.np(1)[1, 3, 3] = nodata
        jim.np(1)[2, 3, 3] = nodata + 1
        jim.np(0)[0, 4, 4] = nodata + 1
        jim.np(1)[0, 4, 4] = nodata
        jim.np(1)[1, 4, 4] = nodata
        jim.np(1)[2, 4, 4] = nodata
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        jim.np(1)[2, 0, 0] = nodata
        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(jim, nodata=nodata, rule='median')
        jim.geometry.reducePlane(nodata=nodata, rule='median')

        stats_reduced = jim.stats.getStats()

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed)'
        assert jim.np()[0, 0] == int(np.median([max, max, min])), \
            'Error in geometry.reducePlane() ' \
            '(rule="mean" did not return median value for all the planes)'
        assert jim.np(1)[0, 0] == int(np.median([3 * max, 2 * max])), \
            'Error in geometry.reducePlane(ref_band) ' \
            '(for rule="median", the used indices of max values are not ' \
            'the ones from the ref_band)'
        assert jim.np(0)[3, 3] == int(np.median([nodata + 1, nodata + 2])), \
            'Error in geometry.reducePlane(nodata) ' \
            '(not ignoring the nodata values for rule="median")'
        assert jim.np(1)[3, 3] == nodata + 1, \
            'Error in geometry.reducePlane(nodata) ' \
            '(not ignoring the nodata values for other bands than band 0 ' \
            'for rule="median")'
        assert jim.np(1)[4, 4] == nodata, \
            'Error in geometry.reducePlane(nodata) ' \
            '(ignoring also the nodata values in other bands than ref_band ' \
            'for rule="median")'
        assert jim.np(0)[2, 2] == nodata, \
            'Error in geometry.reducePlane(nodata) ' \
            '(the returned object not containing nodata on a place where ' \
            'nodata was in all planes for rule="median")'
        assert jim.np(1)[2, 2] == nodata + 1, \
            'Error in geometry.reducePlane(nodata) ' \
            '(the returned object not containing nodata in all bands on a ' \
            'place where nodata was in all planes in the ref_band for ' \
            'rule="median")'
        assert stats_reduced['max'][0] <= stats['max'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the maximum value of returned object is not' \
            ' <= the maximum of the original object)'
        assert stats_reduced['min'][0] >= stats['min'][0], \
            'Error in geometry.reducePlane() ' \
            '(for rule="median", the min value of returned object is not >= ' \
            'the mean of the original object)'

        # Test call with a callback rule
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=3,
                     nband=2, otype='int16', uniform=[min, max])
        jim[2, 0, 0] = max
        jim[1, 0, 0] = max
        jim[0, 0, 0] = min
        jim.np(1)[0, 0, 0] = 3 * max
        jim.np(1)[1, 0, 0] = 2 * max
        jim.np(1)[2, 0, 0] = nodata

        stats = jim.stats.getStats()

        reduced = pj.geometry.reducePlane(
            jim, rule=lambda stacked, plane: stacked * plane)
        jim.geometry.reducePlane(rule=lambda stacked, plane: stacked * plane)

        stats_reduced = jim.stats.getStats(['max', 'min', 'mean'])

        assert jim.properties.isEqual(reduced), \
            'Inconsistency in geometry.reducePlane() ' \
            '(method returns different result than function for a callback ' \
            'rule)'
        assert jim.properties.nrOfPlane() == 1, \
            'Error in geometry.reducePlane() ' \
            '(number of planes not reduced to 1 for a callback rule)'
        assert jim.properties.nrOfBand() == 2, \
            'Error in geometry.reducePlane() ' \
            '(number of bands changed for a callback rule)'
        assert jim.properties.nrOfRow() == jim.properties.nrOfCol() == \
               nr_of_row, \
            'Error in geometry.reducePlane() ' \
            '(number of rows or number of columns changed for a callback rule)'
        assert jim.np()[0, 0] == max * max * min, \
            'Error in geometry.reducePlane() ' \
            '(callback rule not doing what it is supposed to do)'
        assert jim.np(1)[0, 0] == 3 * max * 2 * max * nodata, \
            'Error in geometry.reducePlane() ' \
            '(callback rule not working for all bands)'
        assert stats_reduced['max'][0] >= stats['max'][0], \
            'Error in geometry.reducePlane() ' \
            '(callback rule not doing what it is supposed to do - the max ' \
            'value after using a multiplication rule not >= the max of the ' \
            'original object)'
        assert stats_reduced['min'][0] >= stats['min'][0], \
            'Error in geometry.reducePlane() ' \
            '(callback rule not doing what it is supposed to do - the min ' \
            'value after using a multiplication rule not >= the min of the ' \
            'original object)'
        assert stats_reduced['mean'][0] >= stats['mean'][0], \
            'Error in geometry.reducePlane() ' \
            '(callback rule not doing what it is supposed to do - the mean ' \
            'value after using a multiplication rule not >= the mean of the ' \
            'original object)'

        # Test wrong call with one-plane Jim
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=1, otype='Byte',
                     uniform=[min, max])

        warnings.filterwarnings('error', category=Warning)
        try:
            _ = pj.geometry.reducePlane(jim, 'mean')
            raised = False
        except Warning:
            raised = True

        assert raised, 'Error in raising a warning when performing ' \
                       'reducePlane function on an object with one plane'

        try:
            jim.geometry.reducePlane('mean')
            raised = False
        except Warning:
            raised = True

        assert raised, 'Error in raising a warning when performing ' \
                       'reducePlane method on an object with one plane'

        warnings.resetwarnings()

        # Test call with a non-supported function
        jim = pj.Jim(nrow=nr_of_row, ncol=nr_of_col, nplane=2, otype='Byte',
                     uniform=[min, max])

        try:
            _ = pj.geometry.reducePlane(jim, 'non-supported function')
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, 'Error in raising an error when performing ' \
                       'reducePlane function with a non-supported function'

        try:
            jim.geometry.reducePlane('non-supported function')
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, 'Error in raising an error when performing ' \
                       'reducePlane method with a non-supported function'

        # Test call with max and nodata defined, but no ref_band defined
        try:
            _ = pj.geometry.reducePlane(jim, rule='max', nodata=nodata)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, 'Error in raising an error when performing ' \
                       'reducePlane function with a max/min rule and only ' \
                       'nodata without ref_band defined'

        try:
            jim.geometry.reducePlane(rule='max', nodata=nodata)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, 'Error in raising an error when performing ' \
                       'reducePlane method with a max/min rule and only ' \
                       'nodata without ref_band defined'

        # Test call with rule=callback and nodata specified
        try:
            _ = pj.geometry.reducePlane(
                jim,
                rule=lambda stacked, plane: stacked * plane,
                nodata=nodata)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, 'Error in raising an error when performing ' \
                       'reducePlane function with a callback rule and nodata' \
                       ' defined'

        try:
            jim.geometry.reducePlane(
                rule=lambda stacked, plane: stacked * plane, nodata=nodata)
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, 'Error in raising an error when performing ' \
                       'reducePlane method with a callback rule and nodata' \
                       ' defined'

    @staticmethod
    def test_repeat():
        """Test the repeat() function and method."""

        magnify = 2
        jim = pj.Jim(rasterfn)
        bbox = jim.properties.getBBox()
        ncol = jim.properties.nrOfCol()
        nrow = jim.properties.nrOfRow()
        nplane = jim.properties.nrOfPlane()
        jim.geometry.repeat(magnify,axis=0)
        assert jim.properties.nrOfPlane() == nplane, \
            'Error: number of planes should be identical after repeat 2 row'
        assert jim.properties.nrOfRow() == nrow * magnify, \
            'Error: number of rows should be doubled after repeat 2 row'
        assert jim.properties.nrOfCol() == ncol, \
            'Error: number of cols should be identical after repeat 2 row'
        assert jim.properties.getBBox() == bbox, \
            'Error: bounding box should be identical after repeat 2 row'

        jim.geometry.repeat(magnify,axis=1)
        assert jim.properties.nrOfRow() == nrow * magnify, \
            'Error: number of rows should be doubled after repeat 2 row'
        assert jim.properties.nrOfCol() == ncol * magnify, \
            'Error: number of cols should be doubled after repeat 2 col'
        assert jim.properties.getBBox() == bbox, \
            'Error: bounding box should be identical after repeat 2 col'

        jim = pj.Jim(rasterfn, band2plane=True)
        nplane = jim.properties.nrOfPlane()
        jim.geometry.repeat(magnify,axis=0)
        assert jim.properties.nrOfPlane() == nplane * magnify, \
            'Error: number of planes should be doubled after repeat 2 plane'
        assert jim.properties.nrOfRow() == nrow, \
            'Error: number of rows should be identical after repeat 2 plane'
        assert jim.properties.nrOfCol() == ncol, \
            'Error: number of cols should be identical after repeat 2 plane'
        assert jim.properties.getBBox() == bbox, \
            'Error: bounding box should be identical after repeat 2 plane'

        jim.geometry.repeat(magnify,axis=1)
        assert jim.properties.nrOfPlane() == nplane * magnify, \
            'Error: number of planes should be doubled after repeat 2 plane'
        assert jim.properties.nrOfRow() == nrow * magnify, \
            'Error: number of rows should be doubled after repeat 2 row'
        assert jim.properties.nrOfCol() == ncol, \
            'Error: number of cols should be identical after repeat 2 row'
        assert jim.properties.getBBox() == bbox, \
            'Error: bounding box should be identical after repeat 2 row'

        jim.geometry.repeat(magnify,axis=2)
        assert jim.properties.nrOfPlane() == nplane * magnify, \
            'Error: number of planes should be doubled after repeat 2 plane'
        assert jim.properties.nrOfRow() == nrow * magnify, \
            'Error: number of rows should be doubled after repeat 2 row'
        assert jim.properties.nrOfCol() == ncol * magnify, \
            'Error: number of cols should be doubled after repeat 2 col'
        assert jim.properties.getBBox() == bbox, \
            'Error: bounding box should be identical after repeat 2 col'

class BadGeometryLists(unittest.TestCase):
    """Test functions and methods from geometry module."""

    @staticmethod
    def test_stack():
        """Test the stackBand and stackPlane methods."""
        nrow = ncol = 5
        nplane = nband = 2
        jim1 = pj.Jim(nrow=nrow, ncol=ncol, nband=nband, nplane=nplane,
                      otype='Byte', uniform=10)
        jim2 = pj.Jim(nrow=nrow, ncol=ncol, nband=nband, nplane=nplane,
                      otype='Byte', uniform=10)
        jiml = pj.JimList([jim1, jim2])

        # Test stacking bands of a JimList
        stacked = pj.geometry.stackBand(jiml)
        stacked2 = jiml.geometry.stackBand()

        assert stacked.properties.isEqual(stacked2), \
            'Inconsistency in geometry.stackBand(JimList) ' \
            '(method returns different result than function)'
        assert stacked.properties.nrOfBand() == 2 * nband, \
            'Error in geometry.stackBand(JimList) ' \
            '(number of bands of the returned object not equal to count of ' \
            'bands of all contained Jims)'
        assert (stacked.np(0) == jim1.np(0)).all(), \
            'Error in geometry.stackBand(JimList) ' \
            '(the first band not corresponding to the first band of the ' \
            'first Jim contained)'
        assert (stacked.np(1) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList) ' \
            '(the second band not corresponding to the second band of the ' \
            'first Jim contained)'
        assert (stacked.np(2) == jim2.np(0)).all(), \
            'Error in geometry.stackBand(JimList) ' \
            '(not working for the second Jim contained)'
        assert (stacked.np(3) == jim2.np(1)).all(), \
            'Error in geometry.stackBand(JimList) ' \
            '(not working for other bands than the first one for the second ' \
            'Jim contained)'

        # Test stacking bands of a JimList with band specified
        stacked = pj.geometry.stackBand(jiml, band=1)
        stacked2 = jiml.geometry.stackBand(band=1)

        assert stacked.properties.isEqual(stacked2), \
            'Inconsistency in geometry.stackBand(JimList, band=band) ' \
            '(method returns different result than function)'
        assert stacked.properties.nrOfBand() == 2, \
            'Error in geometry.stackBand(JimList, band=band) ' \
            '(number of bands of the returned object not equal to count of ' \
            'of contained Jims)'
        assert (stacked.np(0) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList, band=band) ' \
            '(the first band not corresponding to the specified band of the ' \
            'first Jim contained)'
        assert (stacked.np(1) == jim2.np(1)).all(), \
            'Error in geometry.stackBand(JimList, band=band) ' \
            '(band x not corresponding to the specified band of the x-th Jim' \
            ' contained for x > 1)'

        # Test stacking bands of two JimLists
        stacked = pj.geometry.stackBand(jiml, jiml)
        stacked2 = jiml.geometry.stackBand(jiml)

        assert stacked.properties.isEqual(stacked2), \
            'Inconsistency in geometry.stackBand(JimList, JimList) ' \
            '(method returns different result than function)'
        assert stacked.properties.nrOfBand() == 2 * 2 * nband, \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(number of bands of the returned object not equal to count of ' \
            'bands of all contained Jims in all JimLists)'
        assert (stacked.np(0) == jim1.np(0)).all(), \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(the first band not corresponding to the first band of the ' \
            'first Jim contained in the first JimList)'
        assert (stacked.np(1) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(the second band not corresponding to the second band of the ' \
            'first Jim contained in the first JimList)'
        assert (stacked.np(2) == jim2.np(0)).all(), \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(not working for the second Jim contained in the first JimList)'
        assert (stacked.np(3) == jim2.np(1)).all(), \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(not working for other bands than the first one for the second ' \
            'Jim contained in the first JimList)'
        assert (stacked.np(4) == jim1.np(0)).all(), \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(not stacking the second JimList after the first one)'
        assert (stacked.np(5) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(not working for other bands than the first one for the second ' \
            'JimList)'
        assert (stacked.np(6) == jim2.np(0)).all(), \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(not working for the second Jim contained in the second JimList)'
        assert (stacked.np(7) == jim2.np(1)).all(), \
            'Error in geometry.stackBand(JimList, JimList) ' \
            '(not working for other bands than the first one for the second ' \
            'Jim contained in the second JimList)'

        # Test stacking bands of two JimLists with band specified
        stacked = pj.geometry.stackBand(jiml, jiml, band=1)
        stacked2 = jiml.geometry.stackBand(jiml, band=1)

        assert stacked.properties.isEqual(stacked2), \
            'Inconsistency in geometry.stackBand(JimList, JimList, ' \
            'band=band) ' \
            '(method returns different result than function)'
        assert stacked.properties.nrOfBand() == 2 * len(jiml), \
            'Error in geometry.stackBand(JimList, JimList, band=band) ' \
            '(number of bands of the returned object not equal to count of ' \
            'of contained Jims)'
        assert (stacked.np(0) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList, JimList, band=band) ' \
            '(the first band not corresponding to the specified band of the ' \
            'first Jim contained)'
        assert (stacked.np(1) == jim2.np(1)).all(), \
            'Error in geometry.stackBand(JimList, JimList, band=band) ' \
            '(band x not corresponding to the specified band of the x-th Jim' \
            ' contained for x > 1)'
        assert (stacked.np(2) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList, JimList, band=band) ' \
            '(not stacking the second JimList after the first one)'
        assert (stacked.np(3) == jim2.np(1)).all(), \
            'Error in geometry.stackBand(JimList, JimList, band=band) ' \
            '(band b+x not corresponding to the specified band of the x-th ' \
            'Jim contained for x > 1 for the second JimList, if b is the sum' \
            ' of bands of the first JimList)'

        # Test stacking bands of a JimList and a Jim
        stacked = pj.geometry.stackBand(jiml, jim1)
        stacked2 = jiml.geometry.stackBand(jim1)

        assert stacked.properties.isEqual(stacked2), \
            'Inconsistency in geometry.stackBand(JimList, Jim) ' \
            '(method returns different result than function)'
        assert stacked.properties.nrOfBand() == 2 * nband + nband, \
            'Error in geometry.stackBand(JimList, Jim) ' \
            '(number of bands of the returned object not equal to count of ' \
            'bands of all contained Jims in all JimLists)'
        assert (stacked.np(0) == jim1.np(0)).all(), \
            'Error in geometry.stackBand(JimList, Jim) ' \
            '(the first band not corresponding to the first band of the ' \
            'first Jim contained in the JimList)'
        assert (stacked.np(1) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList, Jim) ' \
            '(the second band not corresponding to the second band of the ' \
            'first Jim contained in the first JimList)'
        assert (stacked.np(2) == jim2.np(0)).all(), \
            'Error in geometry.stackBand(JimList, Jim) ' \
            '(the second Jim contained in the JimList is not following the ' \
            'first Jim)'
        assert (stacked.np(3) == jim2.np(1)).all(), \
            'Error in geometry.stackBand(JimList, Jim) ' \
            '(not working for other bands than the first one for the second ' \
            'Jim contained in the JimList)'
        assert (stacked.np(4) == jim1.np(0)).all(), \
            'Error in geometry.stackBand(JimList, Jim) ' \
            '(not stacking the Jim after the JimList)'

        # Test stacking bands of a JimList and a Jim with band specified
        stacked = pj.geometry.stackBand(jiml, jim1, band=1)
        stacked2 = jiml.geometry.stackBand(jim1, band=1)

        assert stacked.properties.isEqual(stacked2), \
            'Inconsistency in geometry.stackBand(JimList, Jim, ' \
            'band=band) ' \
            '(method returns different result than function)'
        assert stacked.properties.nrOfBand() == len(jiml) + 1, \
            'Error in geometry.stackBand(JimList, Jim, band=band) ' \
            '(number of bands of the returned object not equal to count of ' \
            'of contained Jims)'
        assert (stacked.np(0) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList, Jim, band=band) ' \
            '(the first band not corresponding to the specified band of the ' \
            'first Jim contained)'
        assert (stacked.np(1) == jim2.np(1)).all(), \
            'Error in geometry.stackBand(JimList, Jim, band=band) ' \
            '(band x not corresponding to the specified band of the x-th Jim' \
            ' contained in the JimList for x > 1)'
        assert (stacked.np(2) == jim1.np(1)).all(), \
            'Error in geometry.stackBand(JimList, Jim, band=band) ' \
            '(not stacking the second JimList after the first one)'

        # Test stacking planes of a JimList
        stacked = pj.geometry.stackPlane(jiml)
        stacked2 = jiml.geometry.stackPlane()

        stacked_first_two_planes = pj.geometry.cropPlane(stacked, [0, 1])
        stacked_last_two_planes = pj.geometry.cropPlane(stacked, [2, 3])

        assert stacked.properties.isEqual(stacked2), \
            'Inconsistency in geometry.stackPlane(JimList) ' \
            '(method returns different result than function)'
        assert stacked.properties.nrOfPlane() == len(jiml) * nplane, \
            'Error in geometry.stackPlane(JimList) ' \
            '(number of planes of the returned object not equal to count of ' \
            'planes of all contained Jims)'
        assert stacked_first_two_planes.properties.isEqual(jim1), \
            'Error in geometry.stackPlane(JimList) ' \
            '(the first planes not corresponding to the planes of the ' \
            'first Jim contained)'
        assert stacked_last_two_planes.properties.isEqual(jim2), \
            'Error in geometry.stackPlane(JimList) ' \
            '(the second Jim contained not stacked after the first Jim)'

        # Test stacking planes of two JimLists
        stacked = pj.geometry.stackPlane(jiml, jiml)
        stacked2 = jiml.geometry.stackPlane(jiml)

        stacked_first_two_planes = pj.geometry.cropPlane(stacked, [0, 1])
        stacked_second_two_planes = pj.geometry.cropPlane(stacked, [2, 3])
        stacked_third_two_planes = pj.geometry.cropPlane(stacked, [4, 5])
        stacked_fourth_two_planes = pj.geometry.cropPlane(stacked, [6, 7])

        assert stacked.properties.isEqual(stacked2), \
            'Inconsistency in geometry.stackPlane(JimList, JimList) ' \
            '(method returns different result than function)'
        assert stacked.properties.nrOfPlane() == 2 * len(jiml) * nplane, \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(number of planes of the returned object not equal to count of ' \
            'planes of all contained Jims in all JimLists)'
        assert stacked_first_two_planes.properties.isEqual(jim1), \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(the first planes not corresponding to the first Jim of the ' \
            'first Jim contained in the first JimList)'
        assert stacked_second_two_planes.properties.isEqual(jim2), \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(not working for the second Jim contained in the first JimList)'
        assert stacked_third_two_planes.properties.isEqual(jim1), \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(not stacking the second JimList after the first one)'
        assert stacked_fourth_two_planes.properties.isEqual(jim2), \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(not working for the second Jim contained in the second JimList)'
        assert (stacked_first_two_planes.np(1) == jim1.np(1)).all(), \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(not working for other bands than the first one)'
        assert (stacked_second_two_planes.np(1) == jim2.np(1)).all(), \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(not working for other bands than the first one for the second ' \
            'Jim contained in the first JimList)'
        assert (stacked_third_two_planes.np(1) == jim1.np(1)).all(), \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(not working for other bands than the first one for the first ' \
            'Jim contained in the second JimList)'
        assert (stacked_fourth_two_planes.np(1) == jim2.np(1)).all(), \
            'Error in geometry.stackPlane(JimList, JimList) ' \
            '(not working for other bands than the first one for the second ' \
            'Jim contained in the second JimList)'

        # Test wrong call with parameter band exceeding the nrOfBand in Jims
        try:
            _ = pj.geometry.stackBand(jiml, jim1,
                                      band=jim1.properties.nrOfBand() + 1)
            raised = False
        except pj.exceptions.JimBandsError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackBand(jim, jim, band) ' \
            'function where the band argument exceeds the number of bands of' \
            ' Jims'

        try:
            jiml.geometry.stackBand(jim1, band=jim1.properties.nrOfBand() + 1)
            raised = False
        except pj.exceptions.JimBandsError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.stackBand(jim, jim, band) ' \
            'method where the band argument exceeds the number of bands of' \
            ' Jims'


class BadGeometryVects(unittest.TestCase):
    """Test functions and methods from geometry module."""

    @staticmethod
    def test_intersect():
        """Test the stack band method."""
        jim = pj.Jim(rasterfn, band=0)
        jimv = pj.JimVect(vectorfn)

        bbox = jim.properties.getBBox()
        new_ulx = (bbox[0] + bbox[2]) / 2.0
        jim_cropped = pj.geometry.crop(jim, ulx=new_ulx, uly=bbox[1],
                                       lrx=bbox[2], lry=bbox[3])

        nr_of_features = jimv.properties.getFeatureCount()

        non_existing_path = pj._get_random_path()

        intersected = pj.geometry.intersect(jimv, jim_cropped,
                                            non_existing_path)
        jimv.geometry.intersect(jim_cropped)

        feature_count_func = intersected.properties.getFeatureCount()
        feature_count_meth = jimv.properties.getFeatureCount()

        assert feature_count_func == feature_count_meth, \
            'Inconsistency in geometry.intersect() ' \
            '(method returns different result than function)'
        assert feature_count_meth < nr_of_features, \
            'Error in geometry.intersect() ' \
            '(not less features on intersected area than on the whole)'

        try:
            _ = pj.geometry.intersect(jimv, jimv, pj._get_random_path())
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.intersect(Jim) ' \
            'function where the argument is not an instance of a Jim object'

        try:
            _ = jimv.geometry.intersect(jimv)
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.intersect(Jim) method ' \
            'where the argument is not an instance of a Jim object'

        intersected.io.close()
        jimv.io.close()

        os.remove(non_existing_path)

    @staticmethod
    def test_convexhull():
        """Test the convexHull() function and method."""
        jimv = pj.JimVect(vectorfn)

        orig_bbox = jimv.properties.getBBox()

        non_existing_path = pj._get_random_path()
        hull = pj.geometry.convexHull(jimv, non_existing_path)
        jimv.geometry.convexHull()

        feature_count_func = hull.properties.getFeatureCount()
        feature_count_meth = jimv.properties.getFeatureCount()

        assert feature_count_func == feature_count_meth, \
            'Inconsistency in geometry.convexHull() ' \
            '(method returns different result than function)'
        assert feature_count_meth == 1, \
            'Error in geometry.convexHull() ' \
            '(getFeatureCount != 1 for the output)'
        assert feature_count_meth == 1, \
            'Error in geometry.convexHull() ' \
            '(getFeatureCount != 1 for the output)'
        assert orig_bbox == jimv.properties.getBBox(), \
            'Error in geometry.convexHull() ' \
            '(BBox of hull is not the same as of the original JimVect)'

        os.remove(non_existing_path)

    @staticmethod
    def test_join():
        """Test the join() function and method."""
        jimv = pj.JimVect(vectorfn)
        jimvid = pj.JimVect(jimv, output='/vsimem/jimvid', newfield='fid',
                            co=['OVERWRITE=YES'])
        jimv.io.close()
        jimr = pj.Jim(rasterfn, band=[0, 1])

        non_existing_path0 = pj._get_random_path()
        non_existing_path1 = pj._get_random_path()
        non_existing_path_joined = pj._get_random_path()

        jimr0 = pj.geometry.cropBand(jimr, 0)
        jimr1 = pj.geometry.cropBand(jimr, 1)

        vect0 = pj.geometry.extract(jimvid, jimr0,
                                    rule='mean',
                                    output=non_existing_path0,
                                    bandname='B0')
        vect1 = pj.geometry.extract(jimvid, jimr1,
                                    rule='mean',
                                    output=non_existing_path1,
                                    bandname='B1')

        vect0_field_names = vect0.properties.getFieldNames()
        vect1_field_names = vect1.properties.getFieldNames()

        joined = pj.geometry.join(vect0, vect1,
                                  output=non_existing_path_joined, fid=['fid'])
        vect0.geometry.join(vect1, output=non_existing_path_joined,
                            fid=['fid'])

        jimvid.io.close()

        feature_count_func = joined.properties.getFeatureCount()
        feature_count_meth = vect0.properties.getFeatureCount()
        bbox_func = joined.properties.getBBox()
        bbox_meth = vect0.properties.getBBox()
        joined_field_names = vect0.properties.getFieldNames()

        assert feature_count_func == feature_count_meth, \
            'Inconsistency in geometry.join() ' \
            '(method returns different result than function)'
        assert bbox_func == bbox_meth, \
            'Inconsistency in geometry.join() ' \
            '(method returns different result than function)'
        assert feature_count_meth == vect1.properties.getFeatureCount(), \
            'Error in geometry.join() (feature count changed)'
        assert len(joined_field_names) > len(vect0_field_names), \
            'Error in geometry.join() ' \
            '(number of fields not raised after join)'
        assert all([i in joined_field_names for i in
                    vect0_field_names + vect1_field_names]), \
            'Error in geometry.join() ' \
            '(not containing all the fields of vectors used for the join)'

        # Test catching wrong calls
        try:
            _ = pj.geometry.join(jimv, jimr, output=non_existing_path_joined)
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.join(JimVect, Jim) ' \
            'function where one of the arguments is not an instance of a ' \
            'JimVect object'

        try:
            _ = jimv.geometry.join(jimr, output=non_existing_path_joined)
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.join(Jim) method ' \
            'where the argument is not an instance of a JimVect object'

        os.remove(non_existing_path0)
        os.remove(non_existing_path1)
        os.remove(non_existing_path_joined)

        # OUTER_FULL
        v1 = pj.JimVect('tests/data/v1.json')
        v2 = pj.JimVect('tests/data/v2.json')
        v3 = pj.JimVect('tests/data/v3.json')

        v1.geometry.join(v2, key=['fid'], method='OUTER_FULL')
        assert len(v1.properties.getFieldNames()) == 3, \
            "Error: join OUTER_FULL field count v1 v2"
        assert v1.properties.getFeatureCount() == 5, \
            "Error: join OUTER_FULL feature count v1 v2"
        v1.geometry.join(v3, key=['fid'], method='OUTER_FULL')
        assert len(v1.properties.getFieldNames()) == 4, \
            "Error: join OUTER_FULL field count v1, v2, and v3"
        assert v1.properties.getFeatureCount() == 6, \
            "Error: join OUTER_FULL feature count v1, v2, and v3"

        # INNER
        v1 = pj.JimVect('tests/data/v1.json')
        v2 = pj.JimVect('tests/data/v2.json')
        v3 = pj.JimVect('tests/data/v3.json')

        v1.geometry.join(v2, key=['fid'], method='INNER')
        assert len(v1.properties.getFieldNames()) == 3, \
            "Error: join OUTER_FULL field count v1 v2"
        assert v1.properties.getFeatureCount() == 3, \
            "Error: join INNER feature count v1 v2"
        v1.geometry.join(v3, key=['fid'], method='INNER')
        assert len(v1.properties.getFieldNames()) == 4, \
            "Error: join OUTER_FULL field count v1, v2, and v3"
        assert v1.properties.getFeatureCount() == 3, \
            "Error: join INNER feature count v1, v2, and v3"

        # OUTER_LEFT
        v1 = pj.JimVect('tests/data/v1.json')
        v2 = pj.JimVect('tests/data/v2.json')
        v3 = pj.JimVect('tests/data/v3.json')

        v1.geometry.join(v2, key=['fid'], method='OUTER_LEFT')
        assert len(v1.properties.getFieldNames()) == 3, \
            "Error: join OUTER_FULL field count v1 v2"
        assert v1.properties.getFeatureCount() == 4, \
            "Error: join OUTER_LEFT feature count v1 v2"
        v1.geometry.join(v3, key=['fid'], method='OUTER_LEFT')
        assert len(v1.properties.getFieldNames()) == 4, \
            "Error: join OUTER_FULL field count v1, v2, and v3"
        assert v1.properties.getFeatureCount() == 4, \
            "Error: join OUTER_LEFT feature count v1, v2, and v3"

        # OUTER_RIGHT
        v1 = pj.JimVect('tests/data/v1.json')
        v2 = pj.JimVect('tests/data/v2.json')
        v3 = pj.JimVect('tests/data/v3.json')

        v1.geometry.join(v2, key=['fid'], method='OUTER_RIGHT')
        assert len(v1.properties.getFieldNames()) == 3, \
            "Error: join OUTER_FULL field count v1 v2"
        assert v1.properties.getFeatureCount() == 4, \
            "Error: join OUTER_RIGHT feature count v1 v2"
        v1.geometry.join(v3, key=['fid'], method='OUTER_RIGHT')
        assert len(v1.properties.getFieldNames()) == 4, \
            "Error: join OUTER_FULL field count v1, v2, and v3"
        assert v1.properties.getFeatureCount() == 4, \
            "Error: join OUTER_RIGHT feature count v1, v2, and v3"

    @staticmethod
    def test_append():
        """Test the append methods."""
        jimv1 = pj.JimVect(nutsfn)
        nfeatures1 = jimv1.properties.getFeatureCount()

        jimv2 = pj.JimVect(nutsfn, ln='milano')
        nfeatures2 = jimv2.properties.getFeatureCount()

        jimv3 = pj.JimVect(nutsfn, ln='lodi')
        nfeatures3 = jimv3.properties.getFeatureCount()

        non_existing_path = pj._get_random_path()
        non_existing_path = os.path.join('/vsimem',
                                         os.path.basename(non_existing_path))
        appended = pj.geometry.append(jimv2, jimv3, non_existing_path,
                                      co=['OVERWRITE=YES'])

        assert nfeatures1 == nfeatures2 + nfeatures3, \
            'Error in opening layers ' \
            '(feature count opening layers)'

        assert jimv1.np(ln=0).all() == appended.np(ln=0).all(), \
            'Error in geometry.append() layer 0' \

        assert jimv1.np(ln=1).all() == appended.np(ln=1).all(), \
            'Error in geometry.append() layer 1' \

        nfeatures23 = nfeatures2 + nfeatures3
        assert appended.properties.getFeatureCount() == nfeatures23, \
            'Error in geometry.append() ' \
            '(feature count append)'

        assert np.append(jimv2.np(ln=0), jimv3.np(ln=0), axis=0).all() == \
            np.append(jimv1.np(ln=0), jimv1.np(ln=1), axis=0).all(), \
            'Error in geometry.append() ' \
            '(append)'

        jimv2.geometry.append(jimv3)

        assert appended.np().all() == jimv2.np().all(), \
            'Error in geometry.append() ' \
            '(function is not equal to method)'

        non_existing_path = pj._get_random_path()
        appended.io.write(non_existing_path)
        os.remove(non_existing_path)

        # Test wrong calls

        jim_raster = pj.Jim(tiles[0])

        try:
            _ = pj.geometry.append(jimv1, jim_raster,
                                   output=pj._get_random_path())
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.append(JimVect, Jim) ' \
            'function where one of arguments is not an instance of a JimVect' \
            ' object'

        try:
            jimv1.geometry.append(jim_raster)
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching a call of geometry.append(JimVect, Jim) ' \
            'method where one of arguments is not an instance of a JimVect ' \
            'object'

        jimv1.io.close()
        jimv2.io.close()
        jimv3.io.close()

    @staticmethod
    def test_sample():
        """Test the sample function."""
        jim0 = pj.Jim(rasterfn, band=[0, 1, 2, 3])
        v = pj.geometry.sample(jim0, random=100, buffer=100, rule=['mean'],
                               output='mem01', oformat='Memory')
        assert v.np().shape == (100, jim0.properties.nrOfBand())
        v.io.close()


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadGeometry),
                  loader.loadTestsFromTestCase(BadGeometryLists),
                  loader.loadTestsFromTestCase(BadGeometryVects)]
    return unittest.TestSuite(suite_list)
