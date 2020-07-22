"""
Test suite for basic methods originating in list() for JimLists.
Author(s): Pieter.Kempeneers@ec.europa.eu,
           Ondrej Pesek,
           Pierre.Soille@ec.europa.eu
Copyright (C) 2018-2020 European Union (Joint Research Centre)

This file is part of pyjeo.

pyjeo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyjeo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyjeo.  If not, see <https://www.gnu.org/licenses/>.
"""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']


class BadBasicMethodLists(unittest.TestCase):
    """Test functions and methods on the root level for JimLists."""

    @staticmethod
    def test_jimlist_creations():
        """Test different ways to create JimList."""
        jim1 = pj.Jim(tiles[0])

        # Test JimList creation from a list of jims
        jiml = pj.JimList([jim1, jim1])

        assert len(jiml) == 2, 'Error in JimList creation from a list of ' \
                               'Jims (wrong number of items afterwards)'
        assert jiml[0].properties.isEqual(jim1) and \
               jiml[1].properties.isEqual(jim1), \
            'Error in JimList creation from a list of Jims ' \
            '(items at positions do not correspond with the used ones)'

        # Test JimList creation from a JimList
        jiml2 = pj.JimList(jiml)

        assert len(jiml2) == 2, 'Error in JimList creation from a JimList ' \
                                '(wrong number of items afterwards)'
        assert jiml2[0].properties.isEqual(jim1) and \
               jiml2[1].properties.isEqual(jim1), \
            'Error in JimList creation from a JimList object ' \
            '(items at positions do not correspond with the used ones)'

        # Test JimList creation from a Jim
        jiml = pj.JimList(jim1)

        assert len(jiml) == 1, 'Error in JimList creation from a Jim object ' \
                               '(wrong number of items afterwards)'
        assert jiml[0].properties.isEqual(jim1), \
            'Error in JimList creation from a Jim object ' \
            '(item at position 0 does not correspond with the used Jim)'

        # Test wrong creation of a JimList (integer used as an argument)
        try:
            jiml = pj.JimList(None)
            raised = False
        except ValueError:
            raised = True

        assert raised, \
            'Error in catching a call of JimList creation with wrong value ' \
            'parsed (integer used instead of allowed argument types)'

    @staticmethod
    def test_list_methods():
        """Test basic methods originating in list() inheritance."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(tiles[1])
        jim3 = pj.Jim(testFile, band=[0, 1])
        jiml1 = pj.JimList([jim1, jim2])
        jiml2 = pj.JimList([jim1, jim2, jim3])

        assert len(jiml1) == 2, 'Error in len(JimList)'

        jiml1.append(jim1)

        assert len(jiml1) == 3, 'Error in len(JimList) or JimList.append()'
        assert len(jiml1) == jiml1._jipjimlist.getSize(), \
            'Error in JimList._set() with argument from_list=True'

        image2 = pj.Jim(jiml1._jipjimlist.getImage(2))

        assert jiml1[2].properties.isEqual(jiml1[0]), \
            'Error in JimList.append()'
        assert not jiml1[2].properties.isEqual(jiml1[1]), \
            'Error in JimList[index]'
        assert jiml1[2].properties.isEqual(image2),\
            'Error in JimList.append()'

        assert jiml1.count(jim1) == 2, 'Error in JimList.count(Jim)'

        jiml1.extend(jiml2)

        assert len(jiml1) == 6, 'Error in JimList.extend()'
        assert len(jiml1) == jiml1._jipjimlist.getSize(), \
            'Error in JimList._set() with argument from_list=True'

        image4 = pj.Jim(jiml1._jipjimlist.getImage(4))

        assert jiml1.count(jim1) == 3, 'Error in JimList.count(Jim)'
        assert jiml1[4].properties.isEqual(image4),\
            'Error in JimList.extend()'

        jiml1.insert(1, jim3)

        assert len(jiml1) == 7, 'Error in JimList.insert()'
        assert len(jiml1) == jiml1._jipjimlist.getSize(), \
            'Error in JimList._set() with argument from_list=True'

        image1 = pj.Jim(jiml1._jipjimlist.getImage(1))

        # TODO: the following line results in segmentation fault...
        # assert jiml1.count(jim3) == 2, 'Error in JimList.insert(Jim)'
        assert jiml1[1].properties.isEqual(image1),\
            'Error in JimList.insert()'

        assert jiml1.index(jim3) == 1, 'Error in JimList.index() or insert()'

        popped = jiml1.pop(jiml1.index(jim3))

        assert len(jiml1) == 6, 'Error in JimList.pop()'
        assert len(jiml1) == jiml1._jipjimlist.getSize(), \
            'Error in JimList._set() with argument from_list=True'
        assert jiml1.count(jim3) == 1, 'Error in JimList.pop(Jim)'
        assert popped.properties.isEqual(jim3), 'Error in JimList.pop()'

        jiml1.remove(jim3)

        assert len(jiml1) == 5, 'Error in JimList.remove()'
        assert len(jiml1) == jiml1._jipjimlist.getSize(), \
            'Error in JimList._set() with argument from_list=True'
        assert jiml1.count(jim3) == 0, 'Error in JimList.remove(Jim)'

        jiml2.reverse()

        image1 = pj.Jim(jiml1._jipjimlist.getImage(1))

        assert jiml2.index(jim1) == 2, 'Error in JimList.reverse()'
        assert jiml1[1].properties.isEqual(image1),\
            'Error in JimList.reverse()'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadBasicMethodLists)]
    return unittest.TestSuite(suite_list)
