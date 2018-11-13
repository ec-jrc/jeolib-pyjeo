"""Run tests for all modules."""

import sys
import unittest

from tests import test_pixops


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    return unittest.TestSuite([test_pixops.load_tests()])


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2).run(load_tests())
    if not result.wasSuccessful():
        sys.exit(1)
