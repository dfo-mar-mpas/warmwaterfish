if(!require(librarian)) install.packages("librarian")
pkgs <- c("targets",
          "sf",
          "arcpullr",
          "leaflet",
          "rgbif",
          "rinat",
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
    name = ns_coast,
    command = {
      st_read("data/NS_coastline_project_Erase1 1.shp")
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
  ),

  tar_target(
    name = rinatdata,
    command = {
      bbox <- st_bbox(simplemar)

      inat_list <- map(
        GBIF_keys$Taxon.Name,
        \(taxon) {
          tryCatch({
            message("Fetching iNaturalist data for: ", taxon)

            result <- get_inat_obs(
              taxon_name = taxon,
              bounds = c(bbox["ymin"], bbox["xmin"], bbox["ymax"], bbox["xmax"]),
              quality = "research",
              geo = TRUE,
              maxresults = 10000
            )

            # get_inat_obs returns a 0-row data frame when no results
            if(nrow(result) == 0) {
              message("  No records found for ", taxon)
              return(NULL)
            }

            message("  Found ", nrow(result), " records for ", taxon)
            return(result)

          }, error = function(e) {
            message("  Error or no records for ", taxon, ": ", e$message)
            return(NULL)
          })
        }
      )


      # Combine all results
      inat_data <- inat_list |>
        bind_rows()

      # Convert to sf object and spatial join
      if(nrow(inat_data) > 0) {
        inat_sf <- inat_data |>
          st_as_sf(
            coords = c("longitude", "latitude"),
            crs = 4326,
            remove = FALSE
          ) |>
          st_join(maritimes, join = st_within)
      } else {
        inat_sf <- NULL
      }

      inat_sf
    }
  )
)


