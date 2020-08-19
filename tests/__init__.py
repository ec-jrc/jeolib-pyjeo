"""Run tests for all modules."""
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

import sys
import unittest

from tests import test_ccops, test_classify, test_demops, test_geometry,\
    test_io, test_jim_basics, test_jimlist_basics, test_jimvect_basics,\
    test_ngbops, test_pixops, test_properties, test_stats


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    # NOTE: In case of raising segmentation faults, maybe there is still the
    #       problem with mialib/jiplib in DEMOps. Try to comment them.

    return unittest.TestSuite([test_ccops.load_tests(),
                               test_classify.load_tests(),
                               test_demops.load_tests(),
                               test_geometry.load_tests(),
                               test_io.load_tests(),
                               test_jim_basics.load_tests(),
                               test_jimlist_basics.load_tests(),
                               test_jimvect_basics.load_tests(),
                               test_ngbops.load_tests(),
                               test_pixops.load_tests(),
                               test_properties.load_tests(),
                               test_stats.load_tests()])


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2).run(load_tests())
    if not result.wasSuccessful():
        sys.exit(1)
