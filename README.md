# Sentinel-2 ingest pipeline

Find the least cloudy Sentinel-2 L2A acquisition within an AOI and export a mosaic as a Cloud Optimized GeoTIFF. Products from one satellite pass are combined to preserve a single observation time for later change analysis.

## Setup and run

Install dependencies in your Python environment:

```powershell
py -m pip install -r requirements.txt
```

A Copernicus Data Space account and S3 credentials are required. Set these environment variables in your terminal or PyCharm Run Configuration. Do not commit credentials.

```powershell
$env:CDSE_S3_ACCESS_KEY = "your-access-key"
$env:CDSE_S3_SECRET_KEY = "your-secret-key"

py s2_pipeline.py --aoi data/inputs/aoi.geojson --start 2026-04-01 --end 2026-05-01 --output data/outputs/example
```

The examples use the Windows Python launcher. On Linux/macOS, use `python` and `export NAME="value"` for environment variables. Rasterio requires GDAL support for JPEG2000 and COG files. Dates are inclusive; the AOI must contain valid Polygon or MultiPolygon geometries with a defined CRS. Running without arguments uses the included AOI, the period above and `data/outputs`.

## Method

- Search the CDSE STAC catalogue and group products by `eopf:datatake_id`. Keep products with complementary footprints, including those for the same tile; remove only identical product IDs.
- Assess SCL on a common 20 m grid in an automatically selected UTM CRS, using nearest-neighbour resampling. Count only pixel centres inside the AOI. Overlaps use the first non-zero SCL value in sorted product-ID order.
- Require 100% observed AOI coverage by default. Rank eligible acquisitions by the combined percentage of clouds (classes 8, 9, 10) and shadows (class 3) among observed pixels. Ties use greater coverage, then acquisition ID.
- Export all 12 L2A spectral bands, applying each asset's STAC scale and offset. Resample 10 m bands by averaging, 20 m bands by nearest neighbour, and 60 m bands by bilinear interpolation. Upsampling does not add spatial detail.
- Mask SCL classes 0, 1, 3, 8, 9 and 10, areas outside the AOI, and pixels missing any required band. Snow and unclassified pixels are retained. Process 512 × 512 windows, one spectral band at a time.

## Outputs

| File | Contents |
| --- | --- |
| `reflectance.tif` | Float32 BOA reflectance COG, NoData -9999. Band order: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B11, B12. B10 is unavailable in L2A. |
| `quality.tif` | SCL and `source_index` bands. Positive source indices refer to the JSON source list; 0 means no source. |
| `reflectance.json` | Source IDs and dates, product/tile counts, quality statistics, grid, production time and code version. |

Cloud/shadow percentages refer to observed pixels; coverage and usable-output percentages refer to all AOI pixels. The detailed assessment report is in `data/cache/acquisition_assessment.json`.

## Cache and options

STAC responses and SCL assessments are cached. Full search parameters, AOI geometry and source IDs are checked before reusing cached data. An identical completed output is skipped on a second run. Raw imagery is not permanently cached, so exporting again can require remote reads.

Use `--refresh-search` to check for new products, `--refresh-assessment` to recalculate quality, and `--overwrite` to replace outputs. For example:

```powershell
py s2_pipeline.py --refresh-search --overwrite
```

Other options: `--cache PATH`, `--crs EPSG:CODE` (projected, metres), and `--min-coverage PERCENT` (default 100). See `--help`. Update `CODE_VERSION` when processing logic changes.

## Errors and limitations

Invalid input, no data, insufficient coverage or no usable output pixels stop the run. Catalogue requests have bounded retries. Failed assessments are recorded; metadata flags selection based on an incomplete candidate comparison. Failed exports require `--overwrite` before retrying in the same folder.

There is no maximum cloud threshold and no multitemporal composite. Cloud boundaries are not buffered, so interpolation can mix neighbouring values near masked pixels. Country-scale processing may require tiling and multiple acquisitions.

## Manual checks

Successful runs covered the supplied regional AOI with multiple products, a small AOI, two disconnected polygons, and an area approximately 50 × 50 km. An identical second run skipped processing. A single-acquisition test reported 23.18% clouds plus shadows and produced masked output. Results were inspected in QGIS.

## Questions

### 1. Besides Sentinel-2 imagery, you now also need to ingest a vector dataset (e.g. a shapefile of mineral deposits) and the PDF publication it is based on. What stays the same in your solution and what changes? How do you keep track of the fact that these three things belong together?

I would keep the common stages: downloading the data, checking its validity, and saving the files and their provenance. For vector data, I would check geometry, attributes and CRS and keep all files belonging to the shapefile, while for the PDF I would save the original document and bibliographic information. I would give the three resources a shared dataset ID and record their relationships in JSON, indicating that the deposit data comes from a particular publication and that the satellite image covers the study area. For a larger collection, I would store metadata and relationships in PostgreSQL while keeping the files in a separate data store.

### 2. This process must run weekly for 200 different AOIs. How do you know what is new and what you already have? What do you record about each run?

I would give each AOI a stable ID and track geometry changes separately from changes to processing settings. The cache in my script reduces repeated work, but weekly updates would require refreshing the catalogue and comparing product IDs and versions with those already recorded. For 200 AOIs, I would store this state in PostgreSQL, using PostGIS for AOI geometries, rather than relying only on JSON files. For each run, I would record the AOI, search period, start and end times, status, code version, new and selected products, output quality and any errors.

### 3. Six months from now someone shows a map built from your output and asks where a specific pixel came from. What do you need to have saved to answer?

First, I would identify the output file used to create the map and any later transformations applied to it. In my current solution, the source_index band in quality.tif points to the source product, whose ID and acquisition date are recorded in JSON. To reproduce the value, I would also need the exact source data version, calibration parameters, output grid, resampling method, masking rules and code version. In an extended solution, I would link this information to the output and pipeline run in PostgreSQL while keeping the rasters as COG files. After resampling, a value may depend on several neighbouring pixels from the source product, so the scene ID alone is not enough to reproduce the calculation.

### 4. What changes when the AOI is a whole country? Where is the bottleneck: memory, network or number of scenes?

I would first check whether a single acquisition can cover the whole country, because my current solution limits the output to one satellite pass. If not, I would allow multiple acquisitions and divide processing, and possibly outputs, into tiles, using a suitable CRS. Block-based processing limits array memory use, but more products increase the number of reads, network transfer and JP2 decoding work. I would identify the bottleneck by measuring elapsed time, RAM, network transfer and disk use rather than assuming that memory or the network is always the limiting factor.

### 5. The best result you can build covers only 70% of the AOI, or is 40% cloudy. Do you return it or fail? What should the user of that data know?

By default, my pipeline requires 100% observed coverage, so a candidate covering only 70% of the AOI would be rejected. The user can explicitly allow a partial result through the `--min-coverage` argument, and the actual coverage is recorded in the metadata. If the best eligible acquisition is 40% cloudy, the current implementation exports it after masking, provided that usable pixels remain, because I have not set a maximum cloud-cover threshold. The user should know the coverage before masking, cloud and shadow percentages, the usable AOI percentage after processing, and whether read failures limited the comparison of candidates.
