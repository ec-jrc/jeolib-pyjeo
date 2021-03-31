# Python program to calculate DEM attributes
# Author(s): Pieter.Kempeneers@ec.europa.eu
#
# Copyright (C) 2021 European Union (Joint Research Centre)
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

from pathlib import Path
import numpy as np
import math
import argparse
import richdem as rd
import pyjeo as pj

parser=argparse.ArgumentParser()

parser.add_argument("-input","--input",help="input path for dem",
                    dest="input",required=True,type=str)
parser.add_argument("-output","--output",help="output path for slope",
                    dest="output",required=True,type=str)
parser.add_argument("-output_dem","--output_dem",help="output path for dem",
                    dest="output_dem",required=False,type=str, default = None)
parser.add_argument("-t_srs","--t_srs",
                    help="target spatial reference system (e.g, epsg:32652)",
                    dest="t_srs",required=False,type=str, default = None)
parser.add_argument("-ulx","--ulx",help="upper left X coordinate in dec deg",
                    dest="ulx",required=False,type=float, default = None)
parser.add_argument("-uly","--uly",help="upper left Y coordinate in dec deg",
                    dest="uly",required=False,type=float, default = None)
parser.add_argument("-lrx","--lrx",help="lower right X coordinate in dec deg",
                    dest="lrx",required=False,type=float, default = None)
parser.add_argument("-lry","--lry",help="lower right Y coordinate in dec deg",
                    dest="lry",required=False,type=float, default = None)
parser.add_argument("-tileindex","--tileindex",help="tileindex to split input",
                    dest="tileindex",required=False,type=int,default=0)
parser.add_argument("-tiletotal","--tiletotal",help="total to split input",
                    dest="tiletotal",required=False,type=int,default=1)
parser.add_argument("-overlap","--overlap",help="overlap to split input",
                    dest="overlap",required=False,type=float,default=0)
parser.add_argument("-attribute","--attribute",
                    help="attribute to calculate [slope_riserun, \
                    slope_percentage, slope_degrees, slope_radians, \
                    aspect, curvature, planform_curvature, \
                    profile_curvature]",
                    dest="attribute", required=False,type=str,
                    default='slope_degrees')

args = parser.parse_args()

demfn = args.input

if (args.ulx is not None and
    args.uly is not None and
    args.lrx is not None and
    args.lry is not None):
    bbox = [args.ulx, args.uly, args.lrx, args.lry]
    jimdem = pj.Jim(Path(demfn), bbox = bbox, tileindex = args.tileindex,
                    tiletotal = args.tiletotal, overlap = args.overlap)
else:
    jimdem = pj.Jim(Path(demfn), tileindex = args.tileindex,
                    tiletotal = args.tiletotal, overlap = args.overlap)

if args.t_srs is not None:
    jimdem.geometry.warp(args.t_srs)

if args.output_dem is not None:
    jimdem.io.write(args.output_dem, co = ['COMPRESS=LZW', 'TILED=YES'])

jimdem.pixops.convert('GDT_Float32')

jimdem[jimdem <= 0] = -9999
dem_richdem  = rd.rdarray(jimdem.np(), no_data=-9999)
dem_richdem.geotransform = jimdem.properties.getGeoTransform()
dem_richdem.projection = jimdem.properties.getProjection()

slope = rd.TerrainAttribute(dem_richdem, attrib=args.attribute)
jimdem.np()[:] = slope
jimdem.properties.setNoDataVals(-9999)

jimdem.io.write(args.output, co = ['COMPRESS=LZW', 'TILED=YES'])
