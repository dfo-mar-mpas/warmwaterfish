if(!require(librarian)) install.packages("librarian")
pkgs <- c("targets",
          "sf",
          "arcpullr",
          "leaflet",
          "rgbif",
          "sf",
          "purrr",
          "dplyr")
librarian::shelf(pkgs)

# Set target options:
tar_option_set(
  packages = c(basename(pkgs))
)

# Define the list of targets:
list(
  tar_target(
    name = maritimes,
    # from https://open.canada.ca/data/en/dataset/f089a3f3-45e9-47de-b1c4-170e9950d8e7
    command = get_spatial_layer("https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/open_data_donnees_ouvertes/eastern_canada_marine_spatial_planning_areas/MapServer/0",
                                where = "NAME_E = 'Scotian Shelf and Bay of Fundy'")),

  tar_target(
    name = simplemar,
    command = {
      d <- 5000
      st_simplify(st_buffer(maritimes,d),dTolerance = d)
    }
  ),

  tar_target(
    name = file_GBIF_keys.csv,
    command = "data/GBIF_keys.csv",
    format = "file"
  ),

  tar_target(
    name = GBIF_keys,
    command = read.csv(file_GBIF_keys.csv)
  ),

  tar_target(
    name = gbifdownload,
    command = {
      taxon_preds <- map2(
        GBIF_keys$Taxon.Level,
        GBIF_keys$GBIF.Key,
        \(x, y) pred_in(paste0(x, "Key"), y)
      )

      occ_download(
        pred_within(st_as_text(simplemar$geoms)),
        pred("hasCoordinate", TRUE),
        pred("hasGeospatialIssue", FALSE),
        do.call(pred_or, taxon_preds),
        format = "SIMPLE_CSV",
        user = Sys.getenv("GBIF_USER"),
        pwd = Sys.getenv("GBIF_PWD"),
        email = Sys.getenv("GBIF_EMAIL")
      )
    }
  ),

  tar_target(
    name = gbifdata,
    command = {
      occ_download_wait(gbifdownload)

      inat <- occ_download_import(occ_download_get(gbifdownload)) |>
        st_as_sf(
          coords = c("decimalLongitude", "decimalLatitude"),
          crs = 4326,
          remove = FALSE
        ) |>
        st_join(maritimes, join = st_within)    }
  ),

  tar_target(
    name = inatdata,
    command = {
      gbifdata |>
        filter(institutionCode == "iNaturalist")}
  )
)


