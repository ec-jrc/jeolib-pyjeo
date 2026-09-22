# Author: Pieter.Kempeneers@ec.europa.eu
#
# Copyright (C) 2025-2026 European Union (Joint Research Centre)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
This program reads as input land-use raster maps and calculates for a given
NUTS3 region per class pixel counts within the edges of urban areas.
It therefore identifies urban areas and their edges, and compute buffers
around the urban areas according to user defined distances. As output, it
writes an enriched land use map that combines the different input rasters
and a csv with the per-class pixel counts for each of the the user defined
distances.
"""

from __future__ import annotations

import argparse
import logging
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyjeo as pj
from osgeo import osr

osr.UseExceptions()

logger = logging.getLogger(__name__)

# "Enriched" cropmap class codes (see module docstring for the full
# description).
CLASS_GRASSLAND = 1
CLASS_HERBACEOUS = 2
CLASS_RESIDENTIAL = 3
CLASS_RESIDENTIAL_BORDER = 4

DEFAULT_DISTANCES = (10, 50, 100, 150)


class CoverageError(RuntimeError):
    """Raised when an input raster does not cover the requested bounding box."""


@dataclass
class Config:
    """Resolved command-line configuration for a single run."""

    outputdir: Path
    tmpdir: Path
    dx: float
    dy: float
    nutsid: str
    distances: list[int]
    cty: str
    gra: str
    her: str
    nuts: str
    urban: str
    verbose: bool

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "Config":
        distances = sorted(args.distance) if args.distance else list(DEFAULT_DISTANCES)
        return cls(
            outputdir=Path(args.outputdir),
            tmpdir=Path(args.tmpdir),
            dx=args.dx,
            dy=args.dy,
            nutsid=args.nutsid,
            distances=distances,
            cty=args.cty,
            gra=args.gra,
            her=args.her,
            nuts=args.nuts,
            urban=args.urban,
            verbose=args.verbose,
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-outputdir", "--outputdir", dest="outputdir", required=True, type=str,
        help="output path",
    )
    parser.add_argument(
        "-tmpdir", "--tmpdir", dest="tmpdir", required=False, type=str,
        default="/scratch/bda", help="Temporary directory",
    )
    parser.add_argument(
        "-dx", "--dx", dest="dx", required=False, type=float, default=10,
        help="spatial resolution in x",
    )
    parser.add_argument(
        "-dy", "--dy", dest="dy", required=False, type=float, default=10,
        help="spatial resolution in y",
    )
    parser.add_argument(
        "-nutsid", "--nutsid", dest="nutsid", required=True, type=str,
        help="load nuts region (e.g., BE335)",
    )
    parser.add_argument(
        "-distance", "--distance", dest="distance", required=False, type=int,
        nargs="+", help="provide a list of (meter) distances to process",
    )
    parser.add_argument(
        "-cty", "--cty", dest="cty", required=False, type=str,
        default="data/crop_type.tif", help="Copernicus HRL crop type map",
    )
    parser.add_argument(
        "-gra", "--gra", dest="gra", required=False, type=str,
        default="data/grass.tif", help="Copernicus HRL grassland",
    )
    parser.add_argument(
        "-her", "--her", dest="her", required=False, type=str,
        default="data/herbaceous.tif", help="Copernicus HRL herbaceous",
    )
    parser.add_argument(
        "-nuts", "--nuts", dest="nuts", required=False, type=str,
        default="data/nuts.shp", help="nuts vector",
    )
    parser.add_argument(
        "-urban", "--urban", dest="urban", required=False, type=str,
        default="data/urban.tif", help="urban map",
    )
    parser.add_argument(
        "--verbose", dest="verbose", action="store_true",
        help="verbose for debugging",
    )
    return parser.parse_args(argv)


def get_extended_bbox(
    bbox: list[float], overlap: float = 5,
    dx: float | None = None, dy: float | None = None,
) -> list[float]:
    """Return ``bbox`` grown on every side by ``overlap`` percent, optionally
    rounded outward to a multiple of the pixel size (``dx``/``dy``)."""
    ulx, uly, lrx, lry = bbox

    extra_x = (overlap / 100) * (lrx - ulx)
    if dx is not None:
        extra_x = math.ceil(extra_x / dx) * dx

    extra_y = (overlap / 100) * (uly - lry)
    if dy is not None:
        extra_y = math.ceil(extra_y / dy) * dy

    return [ulx - extra_x, uly + extra_y, lrx + extra_x, lry - extra_y]


def standardize_wkt(wkt_string: str) -> str:
    """Identify the EPSG code from a raw WKT string and return the official,
    standardized WKT representation for that code.

    Falls back to a normalized (traditional GIS axis order) version of the
    input WKT if no EPSG match is found, and returns an empty string for
    malformed input.
    """
    try:
        srs = osr.SpatialReference()
        # If the string is invalid, this raises a RuntimeError.
        srs.ImportFromWkt(wkt_string)

        # Try to automatically discover the matching EPSG code.
        srs.AutoIdentifyEPSG()
        epsg_code = srs.GetAuthorityCode(None)

        if epsg_code:
            standard_srs = osr.SpatialReference()
            standard_srs.ImportFromEPSG(int(epsg_code))
            standard_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            return standard_srs.ExportToWkt()

        # Fallback for valid custom WKTs without an EPSG registry match.
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        return srs.ExportToWkt()

    except RuntimeError as exc:
        # Catches bad/malformed WKT strings or invalid EPSG code lookups.
        logger.error("Error processing projection: %s", exc)
        return ""


def save_raster(jim, name: str, tmpdir: Path, outputdir: Path, **write_kwargs) -> Path:
    """Write ``jim`` to ``tmpdir`` then move it into ``outputdir``.

    Writing locally first and moving into place avoids partially-written
    output files if ``outputdir`` is a slower/networked filesystem.
    """
    co = write_kwargs.pop("co", ["COMPRESS=LZW", "TILED=YES"])
    tmpfn = tmpdir / name
    outputfn = outputdir / name
    logger.info("write to %s", tmpfn)
    jim.io.write(tmpfn, co=co, **write_kwargs)
    logger.info("moving %s to %s", tmpfn, outputfn)
    shutil.move(tmpfn, outputfn)
    return outputfn


def save_histogram(jim, name: str, tmpdir: Path, outputdir: Path) -> Path:
    """Compute a per-class pixel-count histogram for ``jim`` and save it as CSV."""
    stats = jim.stats.getStats("histogram", nodata=0)
    pdstats = pd.DataFrame(stats)
    pdstats = pdstats.rename(columns={"bin": "class", "histogram": "pixcount"})
    pdstats = pdstats.set_index("class")
    pdstats = pdstats.where(pdstats.pixcount > 0).dropna()

    tmpfn = tmpdir / name
    outputfn = outputdir / name
    pdstats.to_csv(tmpfn)
    logger.info("moving %s to %s", tmpfn, outputfn)
    shutil.move(tmpfn, outputfn)
    return outputfn


def build_cropmap(cfg: Config, bbox: list[float]):
    """Merge the Copernicus HRL crop-type, grassland, and herbaceous layers
    into a single enriched crop-type map."""
    cropmap = pj.Jim(cfg.cty, bbox=bbox, dx=cfg.dx, dy=cfg.dy, align=True)
    if not cropmap:
        raise CoverageError(f"no coverage for cropmap (cty={cfg.cty!r})")

    grassmap = pj.Jim(cfg.gra, bbox=bbox, dx=cfg.dx, dy=cfg.dy, align=True)
    if not grassmap:
        raise CoverageError(f"no coverage for grassmap (gra={cfg.gra!r})")

    herbmap = pj.Jim(cfg.her, bbox=bbox, dx=cfg.dx, dy=cfg.dy, align=True)
    if not herbmap:
        raise CoverageError(f"no coverage for herbmap (her={cfg.her!r})")

    logger.info("bbox cropmap: %s", cropmap.properties.getBBox())

    # Keep values larger or equal than 1000.
    mask = cropmap < 1000
    # Herb pixels will get value CLASS_HERBACEOUS.
    cropmap[(herbmap == 1) & mask] = CLASS_HERBACEOUS
    # Grass supersedes Herb, and will be set to CLASS_GRASSLAND.
    cropmap[(grassmap == 1) & mask] = CLASS_GRASSLAND
    # Recode pixel values equal to 65535 as 0 (no data).
    cropmap[cropmap == 65535] = 0
    return cropmap


def load_urban(cfg: Config, cropmap, bbox: list[float]):
    """Load and, if necessary, reproject the urban land-use map so it lines
    up with ``cropmap``."""
    dx = cropmap.properties.getDeltaX()
    dy = cropmap.properties.getDeltaY()

    urban = pj.Jim(cfg.urban, noread=True)
    urban.properties.setProjection(standardize_wkt(urban.properties.getProjection()))
    cropmap.properties.setProjection(standardize_wkt(cropmap.properties.getProjection()))

    if urban.properties.getProjection() == cropmap.properties.getProjection():
        logger.info("No re-projection of urban map needed")
        # Use the cropmap as the reference layer to align: open and read
        # the urban map at resolution and geographic extent from cropmap.
        urban = pj.Jim(cfg.urban, bbox=bbox, dx=dx, dy=dy)
        if not urban:
            raise CoverageError(f"no coverage for urban map (urban={cfg.urban!r})")
        return urban

    logger.info("Re-projection of urban map needed")
    bbox_extended = get_extended_bbox(bbox, overlap=5, dx=dx, dy=dy)
    logger.debug("bbox_extended: %s", bbox_extended)
    # Opening the raster does not re-project it yet.
    urban = pj.Jim(cfg.urban, bbox=bbox_extended, t_srs="epsg:3035", align=True)
    if not urban:
        raise CoverageError(f"no coverage for re-projected urban map (urban={cfg.urban!r})")
    logger.debug("reproject urban map to LAEA")
    urban.geometry.warp(t_srs="epsg:3035", bbox=bbox, dx=dx, dy=dy)
    return urban


def compute_urban_boundary(urban):
    """Return the ring of residential pixels that are on the edge of an
    urban patch (i.e. urban pixels that are not part of the patch core).

    MSPA uses advanced mathematical morphology operations to classify the
    foreground pixels of a binary image (urban vs. non-urban) into 7
    mutually exclusive structural categories: Core, Islet, Perforation,
    Edge, Loop, Bridge, and Branch. See
    https://pyjeo.readthedocs.io/en/latest/3_reference.html#ccops.segmentBinaryPatterns
    for details on parameters. The input image is ``urban`` with pixels
    characterized as foreground = 2 (urban), background = 1 (non-urban).
    """
    mspa = pj.ccops.segmentBinaryPatterns(urban + 1, 8, size=1, transition=1, internal=1)
    boundary = pj.Jim(urban)
    # Only interested in boundaries of urban areas here (non-core pixels).
    # Cores are represented by values 17 (internal) and 117 (external);
    # boundary is urban area without core, so set these core pixels to 0.
    boundary[mspa == 117] = 0
    boundary[mspa == 17] = 0
    # boundary now has value 1 where we have boundaries and 0 elsewhere.
    return boundary


def process_nuts_region(cfg: Config) -> int:
    """Run the full pestirisk pipeline for a single NUTS region.

    Returns 0 on success, or a non-zero error code on missing coverage.
    """
    cfg.outputdir.mkdir(parents=True, exist_ok=True)
    cfg.tmpdir.mkdir(parents=True, exist_ok=True)

    max_distance = max(cfg.distances)
    attribute_filter = f"NUTS_ID='{cfg.nutsid}'"
    logger.debug("Opening nuts %s", cfg.nuts)
    nuts3 = pj.JimVect(cfg.nuts, attributeFilter=attribute_filter)
    logger.debug("Projection nuts3: %s", nuts3.properties.getProjection())
    nuts3_id = pj.JimVect(nuts3, output="/vsimem/nuts3.sqlite", co=["OVERWRITE=YES"])
    bbox = nuts3_id.properties.getBBox()
    # Provide extra distance for influence in image boundaries.
    bbox[0] -= max_distance
    bbox[1] += max_distance
    bbox[2] += max_distance
    bbox[3] -= max_distance
    nuts3.io.close()

    try:
        cropmap = build_cropmap(cfg, bbox)
    except CoverageError as exc:
        logger.error("%s", exc)
        return 1

    bbox = cropmap.properties.getBBox()

    try:
        urban = load_urban(cfg, cropmap, bbox)
    except CoverageError as exc:
        logger.error("%s", exc)
        return 4

    # Create buffer.
    logger.info("mask urban")
    # We are only interested in residential areas (pixel values == 1).
    urban[urban > 1] = 0
    logger.info("mask urban in cropmap")
    # Code residential areas as CLASS_RESIDENTIAL.
    cropmap[urban] = CLASS_RESIDENTIAL

    logger.info("calculate urban edges")
    boundary = compute_urban_boundary(urban)
    # Code boundary residential areas as CLASS_RESIDENTIAL_BORDER, other
    # pixels retain values. The "enriched" cropmap has now the following
    # pixel values:
    #   0: no valid pixel found in CTY, GRA, nor HER
    #   1: grass based on GRA
    #   2: temporary grassland based on HER
    #   3: residential area
    #   4: residential area border contributing to buffer zone
    #   5 and above: see values from CTY
    cropmap[boundary] = CLASS_RESIDENTIAL_BORDER

    save_raster(cropmap, f"cropmap_{cfg.nutsid}.tif", cfg.tmpdir, cfg.outputdir)

    # Calculate the Euclidean distance in squared pixel values. See
    # https://pyjeo.readthedocs.io/en/latest/3_reference.html#ccops._CCOps.distance2dEuclideanSquared
    # The distance is calculated from the background pixels (value 0) to
    # the foreground pixels (value 1). Here we want the distance from the
    # urban pixels (value 1) to the nearest non-urban pixel, so we negate
    # the image first.
    logger.info("negate")
    urban = urban != 1
    logger.info("calculate squared distance")
    urban.ccops.distance2dEuclideanSquared()
    # urban now contains the square of the distance to the urban areas.
    # Run the statistics for each of the distances to the urban area.
    dy = cropmap.properties.getDeltaY()
    for distance in cfg.distances:
        logger.info("buffer %s", distance)
        distance_px_sq = (distance / dy) ** 2
        logger.info("create mask")
        # Create a binary mask with values 1 near urban areas.
        mask = urban <= distance_px_sq

        logger.info("mask")
        # Only retain crop mask pixels near to urban areas.
        buffered = cropmap[mask]

        # Normalize the WKT projection string to standard traditional GIS
        # order (X, Y).
        buffered.properties.setProjection(
            standardize_wkt(buffered.properties.getProjection())
        )

        # Crop the output to the nuts3 region.
        logger.info("crop the output to the nuts3 region")
        buffered = buffered[nuts3_id]

        # Write buffered cropmap to file only for the last (largest) distance.
        if distance == cfg.distances[-1]:
            save_raster(
                buffered, f"cropmap_buffered_{cfg.nutsid}.tif",
                cfg.tmpdir, cfg.outputdir,
            )

        # Count the pixel statistics in the urban areas for each crop type.
        save_histogram(buffered, f"{cfg.nutsid}_{distance}.csv", cfg.tmpdir, cfg.outputdir)

    # Normalize the WKT projection string to standard traditional GIS order
    # (X, Y).
    cropmap.properties.setProjection(standardize_wkt(cropmap.properties.getProjection()))
    cropmap = cropmap[nuts3_id]
    # Count the pixel statistics in the urban areas for each crop type for
    # the entire crop map, regardless of the distance to urban areas.
    save_histogram(cropmap, f"{cfg.nutsid}.csv", cfg.tmpdir, cfg.outputdir)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cfg = Config.from_args(args)
    logging.basicConfig(
        level=logging.DEBUG if cfg.verbose else logging.INFO,
        format="%(message)s",
    )
    return process_nuts_region(cfg)


if __name__ == "__main__":
    sys.exit(main())
