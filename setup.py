"""Install the pyJEO package."""
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

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from setuptools import find_packages


setup(
    name='pyjeo',
    version='1.1.3',
    author_email='pieter.kempeneers@.ec.europa.eu',
    url='https://github.com/ec-jrc/jeolib-pyjeo',
    description='https://pyjeo.readthedocs.io/',
    license='GPLv3',
    packages=find_packages(exclude=['doc', 'tests']),
    include_package_data=True,
    install_requires=['numpy', 'rioxarray', 'xarray']
)
