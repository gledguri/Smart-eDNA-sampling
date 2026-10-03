# Build the embedded water-only Pacific example from the study prediction extent.
# Usage: Rscript Code/eDNA_pacific_water_v3.R /path/to/ne_10m_ocean.zip
# Source: https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-ocean/
# Natural Earth ocean, 1:10 million, public domain. Keep islands as polygon holes.
suppressPackageStartupMessages(library(sf))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
sf_use_s2(FALSE)
d <- readRDS("Data/pred_data.rds")
# Same extent as w in 02_Data_manip.R; the coordinate units are kilometres.
# Zone 10N is inferred from the manuscript's West Coast location (no RDS CRS tag).
b <- c(range(d$X_utm), range(d$Y_utm)) * 1000
extent <- st_as_sfc(st_bbox(c(xmin=b[1], ymin=b[3], xmax=b[2], ymax=b[4]), crs=32610))
ocean <- st_read(paste0("/vsizip/", normalizePath(args[1]), "/ne_10m_ocean.shp"), quiet=TRUE)
regional <- suppressWarnings(st_crop(st_make_valid(ocean), c(xmin=-127, ymin=38, xmax=-123, ymax=49)))
water <- st_union(st_intersection(st_geometry(st_transform(regional, 32610)), extent))
water <- st_cast(st_collection_extract(st_make_valid(water), "POLYGON", warn=FALSE), "MULTIPOLYGON")
stopifnot(all(st_is_valid(water)), as.numeric(st_area(water)) > 0)
# Segment long offshore extent edges before geographic conversion.
geo <- st_transform(st_segmentize(water, dfMaxLength=20000), 4326)[[1]]
coordinates <- lapply(geo, function(part) lapply(part, function(ring) unname(round(ring, 7))))
payload <- list(
  source=list(name="Natural Earth ocean", scale="1:10 million", version="5.1.1",
    url="https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-ocean/",
    license="public domain", extentSource="Data/pred_data.rds; 02_Data_manip.R",
    projectedCRS="EPSG:32610 (inferred from manuscript location)",
    projectedWaterAreaKm2=as.numeric(st_area(water))/1e6),
  geometry=list(type="MultiPolygon", coordinates=coordinates))
json <- jsonlite::toJSON(payload, auto_unbox=TRUE, digits=7)
writeLines(c("/* Pacific study extent intersected with Natural Earth ocean; rebuilt by eDNA_pacific_water_v3.R. */",
  paste0("const PACIFIC_WATER_DOMAIN_V3 = ",json,";"),
  "if (typeof module !== 'undefined' && module.exports) module.exports = PACIFIC_WATER_DOMAIN_V3;"),
  "Code/eDNA_pacific_water_v3.js")
cat("Water area:",payload$source$projectedWaterAreaKm2,"km²; parts:",length(geo),
    "; rings:",sum(lengths(geo)),"; vertices:",sum(vapply(geo,function(p)sum(vapply(p,nrow,integer(1))),integer(1))),"\n")
