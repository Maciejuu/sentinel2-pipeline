"""Find the least cloudy Sentinel-2 L2A acquisition for an AOI and export a COG mosaic."""

import argparse
import json
import math
import os
from contextlib import ExitStack
from datetime import date, datetime, timezone
from pathlib import Path

# Fixing conflict with PostgreSQL/PostGIS on Windows.
for variable in ["PROJ_LIB", "PROJ_DATA"]:
    if "postgresql" in os.environ.get(variable, "").lower():
        os.environ.pop(variable)

import geopandas as gpd
import numpy as np
import rasterio
from pyproj import CRS
from pystac import ItemCollection
from pystac_client import Client
from pystac_client.exceptions import APIError
from pystac_client.stac_api_io import StacApiIO
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
from rasterio.session import AWSSession
from rasterio.shutil import copy as copy_raster
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window, transform as window_transform
from shapely.geometry import mapping, shape
from shapely.ops import unary_union
from urllib3.util.retry import Retry

CATALOG_URL = "https://stac.dataspace.copernicus.eu/v1/"
COLLECTION_ID = "sentinel-2-l2a"
CODE_VERSION = "0.1.2"
OUTPUT_NODATA = -9999.0

BAND_ASSETS = {
    "B01": "B01_60m", "B02": "B02_10m", "B03": "B03_10m",
    "B04": "B04_10m", "B05": "B05_20m", "B06": "B06_20m",
    "B07": "B07_20m", "B08": "B08_10m", "B8A": "B8A_20m",
    "B09": "B09_60m", "B11": "B11_20m", "B12": "B12_20m",
}

S3_OPTIONS = {
    "AWS_S3_ENDPOINT": "eodata.dataspace.copernicus.eu",
    "AWS_VIRTUAL_HOSTING": False,
    "GDAL_CACHEMAX": 64 * 1024 * 1024,
    "GDAL_HTTP_TIMEOUT": "60",
    "GDAL_HTTP_MAX_RETRY": "2",
    "GDAL_HTTP_RETRY_DELAY": "5",
    "GDAL_NUM_THREADS": "1",
    "OPJ_NUM_THREADS": "1",
}


def parse_args():
    project_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aoi", type=Path, default=project_dir / "data/inputs/aoi.geojson")
    parser.add_argument("--start", default="2026-04-01", help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end", default="2026-05-01", help="End date, inclusive (YYYY-MM-DD)")
    parser.add_argument("--output", type=Path, default=project_dir / "data/outputs")
    parser.add_argument("--cache", type=Path, default=project_dir / "data/cache")
    parser.add_argument("--crs", help="Optional projected CRS in metres, e.g. EPSG:32634")
    parser.add_argument(
        "--min-coverage", type=float, default=100.0,
        help="Minimum required AOI pixel coverage percentage (default: 100.0)",
    )
    parser.add_argument("--refresh-search", action="store_true", help="Ignore STAC search cache")
    parser.add_argument("--refresh-assessment", action="store_true", help="Ignore SCL assessment cache")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    return parser.parse_args()


def write_json(path, data):
    path = Path(path)
    temporary_file = path.with_suffix(".tmp")
    temporary_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
    temporary_file.replace(path)


def read_aoi(aoi_input_file):
    """Load and validate the AOI polygon geometry."""
    if not aoi_input_file.is_file():
        raise FileNotFoundError(f"AOI file not found: {aoi_input_file}")

    gdf = gpd.read_file(aoi_input_file)
    if gdf.empty or gdf.crs is None:
        raise ValueError("The AOI file must be non-empty and have a defined CRS.")
    if gdf.geometry.isna().any() or not gdf.geometry.is_valid.all():
        raise ValueError("The AOI contains missing or invalid geometries.")
    if gdf.geometry.is_empty.any():
        raise ValueError("The AOI contains empty geometries.")
    if not gdf.geom_type.isin(["Polygon", "MultiPolygon"]).all():
        raise ValueError("The AOI must contain only Polygon or MultiPolygon features.")

    search_geometry = gdf.to_crs(epsg=4326)
    aoi_geometry = unary_union(search_geometry.geometry)
    return search_geometry, aoi_geometry


def search_products(query_parameters, cache_dir, refresh_search):
    """Query the CDSE STAC catalogue or load cached items."""
    period_tag = query_parameters["datetime"].replace("/", "_")
    bbox_tag = "_".join(f"{coord:.2f}" for coord in query_parameters["bbox"])
    search_cache_file = cache_dir / f"stac_search_{period_tag}_{bbox_tag}.json"

    items = None
    if search_cache_file.is_file() and not refresh_search:
        try:
            cached_search = json.loads(search_cache_file.read_text(encoding="utf-8"))
            if cached_search.get("query_parameters") == query_parameters:
                items = list(ItemCollection.from_dict(cached_search["items"]))
                print("Loading STAC search results from cache.")
        except (OSError, ValueError, KeyError, TypeError, AttributeError):
            print("Search cache could not be used, querying again...")

    if items is None:
        print("Querying the CDSE STAC catalogue...")
        retry_policy = Retry(
            total=3,
            backoff_factor=5,
            status_forcelist=[429, 502, 503, 504],
            allowed_methods=frozenset(["GET", "POST"]),
            respect_retry_after_header=True,
        )
        stac_io = StacApiIO(max_retries=retry_policy, timeout=60)
        catalog = Client.open(query_parameters["catalog_url"], stac_io=stac_io)
        search = catalog.search(
            collections=[query_parameters["collection"]],
            bbox=query_parameters["bbox"],
            datetime=query_parameters["datetime"],
            limit=100,
        )
        items = list(search.items())
        if items:
            write_json(search_cache_file, {
                "query_parameters": query_parameters,
                "items": ItemCollection(items).to_dict(),
            })

    if not items:
        raise ValueError("No Sentinel-2 products found for the requested AOI and period.")
    print(f"Products found: {len(items)}")
    return items


def group_products(items, aoi_geometry, min_coverage):
    """Group products by acquisition (datatake_id) and check footprint coverage."""
    acq_group = {}
    seen_ids = set()

    for item in items:
        if item.id in seen_ids or item.geometry is None:
            continue
        seen_ids.add(item.id)
        group_id = item.properties.get("eopf:datatake_id")
        if not group_id:
            continue
        if shape(item.geometry).intersects(aoi_geometry):
            acq_group.setdefault(group_id, []).append(item)

    candidate_acquisitions = {}
    excluded_acquisitions = []

    for group_id, group_items in acq_group.items():
        footprint = unary_union([shape(item.geometry) for item in group_items])
        if min_coverage >= 100.0 and not footprint.covers(aoi_geometry):
            excluded_acquisitions.append(group_id)
            print(f"Skipping {group_id}: incomplete catalogue footprint.")
        else:
            candidate_acquisitions[group_id] = group_items

    if not candidate_acquisitions:
        raise ValueError(
            "No acquisition satisfies the footprint coverage requirement. "
            "Try lowering --min-coverage or expanding the date range."
        )
    print(f"Candidate acquisitions for SCL assessment: {len(candidate_acquisitions)}")
    return candidate_acquisitions, excluded_acquisitions


def build_grid(search_geometry, crs_text=None):
    """Define a common 20 m target grid in a metric CRS."""
    resolution = 20
    analysis_crs = CRS.from_user_input(crs_text) if crs_text else search_geometry.estimate_utm_crs()
    if analysis_crs is None or not analysis_crs.is_projected:
        raise ValueError("A projected metric CRS is required for the target grid.")
    if analysis_crs.axis_info[0].unit_name != "metre":
        raise ValueError("The target CRS must use metre units.")

    analysis_aoi = search_geometry.to_crs(analysis_crs)
    analysis_geometry = unary_union(analysis_aoi.geometry)
    min_x, min_y, max_x, max_y = analysis_geometry.bounds

    left = math.floor(min_x / resolution) * resolution
    bottom = math.floor(min_y / resolution) * resolution
    right = math.ceil(max_x / resolution) * resolution
    top = math.ceil(max_y / resolution) * resolution

    grid_width = int(round((right - left) / resolution))
    grid_height = int(round((top - bottom) / resolution))
    grid_transform = from_origin(left, top, resolution, resolution)

    print(f"Target CRS: {analysis_crs} | Grid: {grid_width}x{grid_height} px ({resolution} m)")
    return {
        "crs": analysis_crs,
        "transform": grid_transform,
        "width": grid_width,
        "height": grid_height,
        "aoi": mapping(analysis_geometry),
        "resolution": resolution,
    }


def assess_acquisitions(candidate_acquisitions, grid, s3_session, cache_dir, refresh_assessment):
    """Evaluate AOI cloud and shadow percentages using windowed SCL_20m reads."""
    block_size = 512 # Limit array memory use.
    acquisition_results = []
    failed_acquisitions = []

    grid_info = {
        "crs": grid["crs"].to_string(),
        "transform": list(grid["transform"])[:6],
        "width": grid["width"],
        "height": grid["height"],
        "aoi": json.loads(json.dumps(grid["aoi"])),
        "code_version": CODE_VERSION,
    }

    with rasterio.Env(session=s3_session, **S3_OPTIONS):
        for group_id, group_items in candidate_acquisitions.items():
            ordered_items = sorted(group_items, key=lambda item: item.id)
            cache_file = cache_dir / f"scl_{group_id}.json"

            if cache_file.is_file() and not refresh_assessment:
                try:
                    result = json.loads(cache_file.read_text(encoding="utf-8"))
                    if (
                        result.get("grid_info") == grid_info
                        and result.get("source_scene_ids") == [item.id for item in ordered_items]
                    ):
                        acquisition_results.append(result)
                        print(f"Loaded cached SCL assessment: {group_id}")
                        continue
                except (OSError, ValueError, KeyError):
                    print(f"Invalid cache for {group_id}, recalculating...")

            print(f"Assessing {group_id} ({len(ordered_items)} products)...", flush=True)
            try:
                with ExitStack() as stack:
                    scl_sources = []
                    for item in ordered_items:
                        asset = item.assets.get("SCL_20m")
                        if asset is None:
                            raise ValueError(f"Missing SCL_20m asset in {item.id}.")
                        src = stack.enter_context(rasterio.open(asset.href))
                        source_crs = src.crs or asset.extra_fields.get("proj:code")
                        vrt = stack.enter_context(WarpedVRT(
                            src,
                            src_crs=source_crs,
                            src_nodata=0,
                            crs=grid["crs"],
                            transform=grid["transform"],
                            width=grid["width"],
                            height=grid["height"],
                            nodata=0,
                            resampling=Resampling.nearest,
                        ))
                        scl_sources.append(vrt)

                    aoi_count = 0
                    observed_count = 0
                    cloud_count = 0
                    shadow_count = 0
                    defective_count = 0
                    usable_count = 0

                    for row in range(0, grid["height"], block_size):
                        for col in range(0, grid["width"], block_size):
                            height = min(block_size, grid["height"] - row)
                            width = min(block_size, grid["width"] - col)
                            window = Window(col, row, width, height)

                            inside_aoi = geometry_mask(
                                [grid["aoi"]],
                                out_shape=(height, width),
                                transform=window_transform(window, grid["transform"]),
                                invert=True,
                                all_touched=False,
                            )
                            if not inside_aoi.any():
                                continue

                            mosaic = np.zeros((height, width), dtype="uint8")
                            for vrt in scl_sources:
                                if not np.any(inside_aoi & (mosaic == 0)):
                                    break
                                data = vrt.read(1, window=window)
                                take = inside_aoi & (mosaic == 0) & (data != 0)
                                mosaic[take] = data[take]

                            observed = inside_aoi & (mosaic != 0)
                            clouds = observed & np.isin(mosaic, [8, 9, 10])
                            shadows = observed & (mosaic == 3)
                            defective = observed & (mosaic == 1)
                            usable = observed & ~clouds & ~shadows & ~defective

                            aoi_count += int(np.count_nonzero(inside_aoi))
                            observed_count += int(np.count_nonzero(observed))
                            cloud_count += int(np.count_nonzero(clouds))
                            shadow_count += int(np.count_nonzero(shadows))
                            defective_count += int(np.count_nonzero(defective))
                            usable_count += int(np.count_nonzero(usable))

                if aoi_count == 0:
                    raise ValueError("AOI polygon covers 0 pixel centres on the 20 m grid.")

                coverage_percent = 100.0 * observed_count / aoi_count
                affected_percent = (
                    100.0 * (cloud_count + shadow_count) / observed_count
                    if observed_count > 0 else None
                )

                result = {
                    "acquisition_id": group_id,
                    "grid_info": grid_info,
                    "source_scene_ids": [item.id for item in ordered_items],
                    "source_dates": [item.datetime.isoformat() for item in ordered_items],
                    "aoi_pixels": aoi_count,
                    "observed_pixels": observed_count,
                    "cloud_pixels": cloud_count,
                    "shadow_pixels": shadow_count,
                    "defective_pixels": defective_count,
                    "usable_pixels": usable_count,
                    "coverage_percent": coverage_percent,
                    "cloud_shadow_percent_observed": affected_percent,
                    "usable_percent_aoi": 100.0 * usable_count / aoi_count,
                    "full_coverage": observed_count == aoi_count,
                }
                acquisition_results.append(result)
                write_json(cache_file, result)

                cs_str = "N/A" if affected_percent is None else f"{affected_percent:.2f}%"
                print(
                    f"  Coverage: {coverage_percent:.2f}% | "
                    f"Clouds+shadows: {cs_str} | "
                    f"Usable AOI: {result['usable_percent_aoi']:.2f}%",
                    flush=True,
                )

            except (rasterio.errors.RasterioError, ValueError, KeyError) as error:
                failed_acquisitions.append({"acquisition_id": group_id, "error": str(error)})
                print(f"  Warning: assessment failed for {group_id}: {error}", flush=True)

    return acquisition_results, failed_acquisitions


def export_bands(best_items, best_result, grid, s3_session, output_dir, query_parameters,
                 aoi_geometry, request_params, min_coverage=100.0, failed_acquisition_count=0):
    """Export the 12-band reflectance COG, 2-band quality COG, and JSON metadata."""
    best_items = sorted(best_items, key=lambda item: item.id)
    best_group_id = best_result["acquisition_id"]

    image_file = output_dir / "reflectance.tif"
    quality_file = output_dir / "quality.tif"
    metadata_file = output_dir / "reflectance.json"
    working_image = output_dir / "reflectance_work.tif"
    working_quality = output_dir / "quality_work.tif"

    base_profile = {
        "driver": "GTiff",
        "width": grid["width"],
        "height": grid["height"],
        "crs": grid["crs"],
        "transform": grid["transform"],
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "compress": "deflate",
        "BIGTIFF": "IF_SAFER",
    }
    quality_profile = {**base_profile, "count": 2, "dtype": "uint16", "nodata": 0}
    image_profile = {
        **base_profile,
        "count": len(BAND_ASSETS),
        "dtype": "float32",
        "nodata": OUTPUT_NODATA,
        "predictor": 3,
    }

    try:
        with rasterio.Env(session=s3_session, **S3_OPTIONS):
            print("Building SCL mosaic and per-pixel source index...", flush=True)
            with ExitStack() as stack:
                scl_sources = []
                for item in best_items:
                    asset = item.assets["SCL_20m"]
                    src = stack.enter_context(rasterio.open(asset.href))
                    source_crs = src.crs or asset.extra_fields.get("proj:code")
                    vrt = stack.enter_context(WarpedVRT(
                        src, src_crs=source_crs, src_nodata=0,
                        crs=grid["crs"], transform=grid["transform"],
                        width=grid["width"], height=grid["height"], nodata=0,
                        resampling=Resampling.nearest,
                    ))
                    scl_sources.append(vrt)

                with rasterio.open(working_quality, "w", **quality_profile) as dst:
                    dst.set_band_description(1, "SCL")
                    dst.set_band_description(2, "source_index")
                    for _, window in dst.block_windows(1):
                        height, width = int(window.height), int(window.width)
                        inside_aoi = geometry_mask(
                            [grid["aoi"]], out_shape=(height, width),
                            transform=window_transform(window, grid["transform"]),
                            invert=True, all_touched=False,
                        )
                        scl = np.zeros((height, width), dtype="uint16")
                        source_map = np.zeros((height, width), dtype="uint16")
                        if inside_aoi.any():
                            for source_index, vrt in enumerate(scl_sources, start=1):
                                if not np.any(inside_aoi & (scl == 0)):
                                    break
                                data = vrt.read(1, window=window)
                                take = inside_aoi & (scl == 0) & (data != 0)
                                scl[take] = data[take]
                                source_map[take] = source_index
                        dst.write(scl, 1, window=window)
                        dst.write(source_map, 2, window=window)

            missing_by_band = {}
            valid_output_count = 0
            source_pixel_counts = [0] * len(best_items)

            with rasterio.open(working_quality) as quality_src:
                with rasterio.open(working_image, "w", **image_profile) as dst:
                    for band_index, (band_name, asset_name) in enumerate(BAND_ASSETS.items(), start=1):
                        print(f"Exporting {band_name} ({band_index}/{len(BAND_ASSETS)})...", flush=True)
                        missing_count = 0
                        dst.set_band_description(band_index, band_name)

                        with ExitStack() as stack:
                            band_sources = []
                            for item in best_items:
                                asset = item.assets[asset_name]
                                scale = asset.extra_fields.get("raster:scale")
                                offset = asset.extra_fields.get("raster:offset")
                                if scale is None or offset is None:
                                    raise ValueError(
                                        f"Missing calibration metadata: {item.id}, {asset_name}."
                                    )
                                scale = float(scale)
                                offset = float(offset)
                                source_nodata = asset.extra_fields.get("nodata", 0)
                                gsd = asset.extra_fields.get("gsd", 20)

                                if gsd < grid["resolution"]:
                                    resampling = Resampling.average
                                elif gsd > grid["resolution"]:
                                    resampling = Resampling.bilinear
                                else:
                                    resampling = Resampling.nearest

                                src = stack.enter_context(rasterio.open(asset.href))
                                source_crs = src.crs or asset.extra_fields.get("proj:code")
                                vrt = stack.enter_context(WarpedVRT(
                                    src, src_crs=source_crs, src_nodata=source_nodata,
                                    crs=grid["crs"], transform=grid["transform"],
                                    width=grid["width"], height=grid["height"],
                                    dtype="float32", nodata=np.nan,
                                    resampling=resampling,
                                ))
                                band_sources.append((vrt, scale, offset))

                            for _, window in dst.block_windows(1):
                                scl = quality_src.read(1, window=window)
                                source_map = quality_src.read(2, window=window)
                                good = (source_map != 0) & ~np.isin(scl, [0, 1, 3, 8, 9, 10])
                                block = np.full(scl.shape, OUTPUT_NODATA, dtype="float32")

                                for source_index, (vrt, scale, offset) in enumerate(band_sources, start=1):
                                    selected = good & (source_map == source_index)
                                    if not selected.any():
                                        continue
                                    raw = vrt.read(1, window=window)
                                    valid = selected & np.isfinite(raw)
                                    block[valid] = raw[valid] * scale + offset
                                    missing_count += int(np.count_nonzero(selected & ~np.isfinite(raw)))

                                dst.write(block, band_index, window=window)
                                if band_index == 1 and missing_count == 0:
                                    valid_first = block != OUTPUT_NODATA
                                    valid_output_count += int(np.count_nonzero(valid_first))
                                    for idx in range(1, len(best_items) + 1):
                                        source_pixel_counts[idx - 1] += int(
                                            np.count_nonzero(valid_first & (source_map == idx))
                                        )
                        missing_by_band[band_name] = missing_count

            if any(missing_by_band.values()):
                valid_output_count = 0
                source_pixel_counts = [0] * len(best_items)
                with rasterio.open(working_image, "r+") as dst, rasterio.open(working_quality) as quality_src:
                    for _, window in dst.block_windows(1):
                        data = dst.read(window=window)
                        valid = np.all(np.isfinite(data) & (data != OUTPUT_NODATA), axis=0)
                        data[:, ~valid] = OUTPUT_NODATA
                        dst.write(data, window=window)
                        valid_output_count += int(np.count_nonzero(valid))
                        source_map = quality_src.read(2, window=window)
                        for idx in range(1, len(best_items) + 1):
                            source_pixel_counts[idx - 1] += int(
                                np.count_nonzero(valid & (source_map == idx))
                            )

            if valid_output_count == 0:
                raise ValueError("No valid cloud-free pixels remain in the output image.")

            print("Converting outputs to Cloud Optimized GeoTIFF...", flush=True)
            metadata_file.unlink(missing_ok=True)
            copy_raster(
                working_image, image_file, driver="COG",
                COMPRESS="DEFLATE", PREDICTOR="FLOATING_POINT",
                BLOCKSIZE=512, BIGTIFF="IF_SAFER", OVERVIEW_RESAMPLING="AVERAGE",
            )
            copy_raster(
                working_quality, quality_file, driver="COG",
                COMPRESS="DEFLATE", BLOCKSIZE=512, BIGTIFF="IF_SAFER",
                OVERVIEW_RESAMPLING="NEAREST",
            )
    finally:
        working_image.unlink(missing_ok=True)
        working_quality.unlink(missing_ok=True)

    sources = [
        {
            "source_index": idx,
            "scene_id": item.id,
            "date": item.datetime.isoformat(),
            "tile": item.properties.get("grid:code"),
            "output_pixels": source_pixel_counts[idx - 1],
        }
        for idx, item in enumerate(best_items, start=1)
    ]

    observed = best_result["observed_pixels"]
    metadata = {
        "produced_at": datetime.now(timezone.utc).isoformat(),
        "code_version": CODE_VERSION,
        "request_params": request_params,
        "acquisition_id": best_group_id,
        "image": image_file.name,
        "quality": quality_file.name,
        "source_product_count": len(sources),
        "contributing_product_count": sum(count > 0 for count in source_pixel_counts),
        "source_tile_count": len({source["tile"] for source in sources}),
        "sources": sources,
        "min_coverage_percent": min_coverage,
        "selection_complete": failed_acquisition_count == 0,
        "failed_acquisition_count": failed_acquisition_count,
        "search_period": query_parameters["datetime"],
        "aoi": mapping(aoi_geometry),
        "quality_percent": {
            "coverage": round(best_result["coverage_percent"], 2),
            "clouds_of_observed": round(100.0 * best_result["cloud_pixels"] / observed, 2),
            "shadows_of_observed": round(100.0 * best_result["shadow_pixels"] / observed, 2),
            "clouds_and_shadows_of_observed": round(best_result["cloud_shadow_percent_observed"], 2),
            "usable_output_of_aoi": round(100.0 * valid_output_count / best_result["aoi_pixels"], 2),
        },
        "raster": {
            "crs": grid["crs"].to_string(),
            "resolution_m": grid["resolution"],
            "width": grid["width"],
            "height": grid["height"],
            "transform": list(grid["transform"])[:6],
            "bands": list(BAND_ASSETS.keys()),
            "dtype": "float32",
            "nodata": OUTPUT_NODATA,
            "values": "BOA surface reflectance",
            "units": "dimensionless; output scale 1, offset 0",
        },
    }
    write_json(metadata_file, metadata)
    print(f"Saved reflectance COG: {image_file}")
    print(f"Saved quality COG: {quality_file}")
    print(f"Saved metadata: {metadata_file}")


def main():
    args = parse_args()
    if date.fromisoformat(args.start) > date.fromisoformat(args.end):
        raise ValueError("Start date must not be later than end date.")
    if not (0.0 < args.min_coverage <= 100.0):
        raise ValueError("--min-coverage must be in the range (0, 100].")

    output_dir = args.output.resolve()
    cache_dir = args.cache.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    search_geometry, aoi_geometry = read_aoi(args.aoi.resolve())
    grid = build_grid(search_geometry, args.crs)

    query_parameters = {
        "catalog_url": CATALOG_URL,
        "collection": COLLECTION_ID,
        "bbox": search_geometry.total_bounds.tolist(),
        "datetime": f"{args.start}/{args.end}",
    }
    request_params = {
        "aoi": json.loads(json.dumps(mapping(aoi_geometry))),
        "aoi_bbox": search_geometry.total_bounds.tolist(),
        "datetime": query_parameters["datetime"],
        "crs": grid["crs"].to_string(),
        "min_coverage": args.min_coverage,
        "code_version": CODE_VERSION,
    }

    metadata_file = output_dir / "reflectance.json"
    image_file = output_dir / "reflectance.tif"
    quality_file = output_dir / "quality.tif"

    output_files = [metadata_file, image_file, quality_file]
    if any(path.exists() for path in output_files) and not args.overwrite:
        if (
            all(path.is_file() for path in output_files)
            and not args.refresh_search and not args.refresh_assessment
        ):
            try:
                existing_meta = json.loads(metadata_file.read_text(encoding="utf-8"))
                if (
                    existing_meta.get("request_params") == request_params
                    and existing_meta.get("selection_complete") is True
                ):
                    print(f"Matching output already exists in {output_dir}. Skipping run.")
                    return
            except (OSError, ValueError):
                pass
        raise FileExistsError(
            f"Output files already exist in {output_dir}. Use --overwrite or a different --output path."
        )

    items = search_products(query_parameters, cache_dir, args.refresh_search)
    candidate_acquisitions, excluded_acquisitions = group_products(
        items, aoi_geometry, args.min_coverage
    )

    access_key = os.environ.get("CDSE_S3_ACCESS_KEY")
    secret_key = os.environ.get("CDSE_S3_SECRET_KEY")
    if not access_key or not secret_key:
        raise ValueError("Missing CDSE_S3_ACCESS_KEY or CDSE_S3_SECRET_KEY environment variables.")

    s3_session = AWSSession(
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        region_name="default",
    )

    acquisition_results, failed_acquisitions = assess_acquisitions(
        candidate_acquisitions, grid, s3_session, cache_dir, args.refresh_assessment
    )

    eligible_results = [
        r for r in acquisition_results
        if r["coverage_percent"] >= args.min_coverage and r["cloud_shadow_percent_observed"] is not None
    ]
    best_result = None
    if eligible_results:
        best_result = min(
            eligible_results,
            key=lambda r: (r["cloud_shadow_percent_observed"], -r["coverage_percent"], r["acquisition_id"]),
        )

    write_json(cache_dir / "acquisition_assessment.json", {
        "produced_at": datetime.now(timezone.utc).isoformat(),
        "code_version": CODE_VERSION,
        "request_params": request_params,
        "results": acquisition_results,
        "excluded_catalogue_coverage": excluded_acquisitions,
        "failed_acquisitions": failed_acquisitions,
        "best_acquisition_id": best_result["acquisition_id"] if best_result else None,
    })

    if failed_acquisitions:
        print(
            f"Warning: {len(failed_acquisitions)} acquisition(s) failed during SCL read "
            "(details saved in cache/acquisition_assessment.json)."
        )

    if best_result is None:
        raise ValueError(
            f"No acquisition met the >= {args.min_coverage:.1f}% pixel coverage requirement. "
            "Try lowering --min-coverage or changing the search period."
        )

    best_group_id = best_result["acquisition_id"]
    print(
        f"Selected acquisition: {best_group_id} "
        f"(coverage: {best_result['coverage_percent']:.2f}%, "
        f"clouds+shadows: {best_result['cloud_shadow_percent_observed']:.2f}%)"
    )

    export_bands(
        candidate_acquisitions[best_group_id],
        best_result,
        grid,
        s3_session,
        output_dir,
        query_parameters,
        aoi_geometry,
        request_params,
        min_coverage=args.min_coverage,
        failed_acquisition_count=len(failed_acquisitions),
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, APIError, rasterio.errors.RasterioError) as error:
        raise SystemExit(f"Pipeline error: {error}") from None