import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


from setuptools import find_packages


setup(
    name='pyjeo',
    version='0.5.0',
    author_email='ondrej.pesek@ec.europa.eu',
    url='https://cidportal.jrc.ec.europa.eu/apps/gitlab/JIPlib/pyJEO',
    packages=find_packages(exclude=['doc', 'tests']),
    include_package_data=True)


