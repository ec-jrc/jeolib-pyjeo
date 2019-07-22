"""Run tests for all modules."""

import sys
import unittest

from tests import test_ccops, test_demops, test_geometry, test_io, \
    test_jim_basics, test_jimlist_basics, test_jimvect_basics, test_ngbops, \
    test_pixops, test_properties, test_stats


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    # NOTE: In case of raising segmentation faults, maybe there is still the
    #       problem with mialib/jiplib in DEMOps. Try to comment them.

    return unittest.TestSuite([test_ccops.load_tests(),
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
