"""Install the pyJEO package."""

import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


from setuptools import find_packages


setup(
    name='pyjeo',
    version='0.5.0',
    author_email='ondrej.pesek@ext.ec.europa.eu',
    url='https://cidportal.jrc.ec.europa.eu/apps/gitlab/JIPlib/pyJEO',
    description='https://cidportal.jrc.ec.europa.eu/'
                'services/processing/jiplibhelp/',
    license='EUPL',
    packages=find_packages(exclude=['doc', 'tests']),
    include_package_data=True,
    install_requires='numpy'
)
