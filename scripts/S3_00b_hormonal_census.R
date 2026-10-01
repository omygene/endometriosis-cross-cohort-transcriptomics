# ============================================================================
# S3_00b_hormonal_census.R -- SUPERSEDED 2026-09-19 by S3_00c (parser v3).
# Kept for provenance only. DO NOT RUN.
# Measured finding that retired this script: GSE179640_filtered.rds carries
# only orig.ident + QC metrics (no treatment metadata, no 'sample' column);
# the hormonal source is the official GEO series matrix, extracted by
# S3_00c_geo_hormonal_census.R. See S3_design.md sections 3, 6 (item 4) and 7.5.
# ============================================================================
set.seed(42)
stop("S3_00b is SUPERSEDED -- do not run. Use scripts/S3_00c_geo_hormonal_census.R ",
     "(see S3_design.md). This stop is intentional (R9: no silent dead code).")
