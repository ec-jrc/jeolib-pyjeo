# HTCondor script to submit dem.py
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

Universe = docker
Executable = /usr/bin/python3
Arguments = dem.py -input /eos/jeodpp/data/base/Elevation/GLOBAL/AW3D30/VER2-1/Data/VRT/aw3d30_dsm_10deg_relpath.vrt -output_dem /eos/jeodpp/data/projects/BDA/training/slope_dem_$(Step)_$(tiletotal)_dem.tif -output /eos/jeodpp/data/projects/BDA/training/slope_dem_$(Step)_$(tiletotal).tif -tileindex $(Step) -tiletotal $(tiletotal) -ulx 5.0 -uly 45.5 -lrx 9.0 -lry 41.5 -t_srs epsg:3035
transfer_input_files = dem.py
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
Docker_image = jeoreg.cidsn.jrc.it:5000/jeodpp-htcondor/base_gdal_py3_deb10_pyjeo:0.5.31
request_cpus = 1
batch_name = dem_training
Priority = 1001
request_memory = 8GB
Requirements = EOSisRunning
Output = /eos/jeodpp/htcondor/processing_logs/BDA/kempepi/training/dem_$(ClusterId)_$(Step).out
Error = /eos/jeodpp/htcondor/processing_logs/BDA/kempepi/training/dem_$(ClusterId)_$(Step).err
Log = /eos/jeodpp/htcondor/processing_logs/BDA/kempepi/training/dem_$(ClusterId)_$(Step).log
queue $(tiletotal)
