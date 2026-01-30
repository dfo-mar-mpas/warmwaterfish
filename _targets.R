if(!require(librarian)) install.packages("librarian")
pkgs <- c("targets",
          "sf",
          "arcpullr",
          "leaflet",
          "rgbif",
          "rinat",
          "sf",
          "purrr",
          "rnaturalearth",
          "worrms",
          "qs",
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
      st_read("data/NS_coastline_project_Erase1 1.shp") |>
        st_transform(st_crs(maritimes))
    }
  ),

  tar_target(
    name = nb_coast,
    command = {
      ne_states(country = "Canada", returnclass = "sf") |>
        filter(name == "New Brunswick")
    }
  ),

  tar_target(
    name = coastalmar,
    command = {
      coast <- st_union(c(ns_coast$geometry,
                          nb_coast$geometry)) |>
        st_make_valid()

      st_intersection(maritimes,
                      coast |>
                        st_buffer(10000))
    }
  ),

  tar_target(
    name = simplecoastalmar,
    command = {
      d <- 5000
      st_simplify(st_buffer(coastalmar,d),dTolerance = d)
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
    name = gbifdownloadall,
    command = {
      occ_download(
        pred_within(st_as_text(simplecoastalmar$geoms)),
        pred("hasCoordinate", TRUE),
        pred("hasGeospatialIssue", FALSE),
        format = "SIMPLE_CSV",
        user = Sys.getenv("GBIF_USER"),
        pwd = Sys.getenv("GBIF_PWD"),
        email = Sys.getenv("GBIF_EMAIL")
      )
    }
  ),

  tar_target(
    name = gbifdataall,
    command = {
      occ_download_wait(gbifdownloadall)

      gbif <- occ_download_import(occ_download_get(gbifdownloadall)) |>
        st_as_sf(
          coords = c("decimalLongitude", "decimalLatitude"),
          crs = 4326,
          remove = FALSE
        ) |>
        st_join(coastalmar, join = st_within) |>
        select(where(~n_distinct(., na.rm = TRUE) > 1)) |> # remove useless columns
        mutate(
          # Determine which name to use based on taxonomic rank and availability
          name_to_check = case_when(
            # If we have species-level ID
            taxonRank == "SPECIES" & !is.na(species) ~ species,
            taxonRank == "SPECIES" & is.na(species) & !is.na(scientificName) ~ scientificName,

            # If identified to genus
            taxonRank == "GENUS" & !is.na(genus) ~ genus,

            # If identified to family
            taxonRank == "FAMILY" & !is.na(family) ~ family,

            # Fallback: use best available name
            !is.na(species) ~ species,
            !is.na(genus) ~ genus,
            !is.na(family) ~ family,
            !is.na(scientificName) ~ scientificName,

            TRUE ~ NA_character_
          )
        )



      marinekey <- gbif |>
        pull(name_to_check) |>
        unique() |>
        sort() |>
        lapply(function(taxon_name){
          if (is.na(taxon_name) || taxon_name == "") {
            return(data.frame(
              is_marine = NA,
              is_brackish = NA,
              is_freshwater = NA,
              is_terrestrial = NA,
              worms_rank = NA,
              aphia_id = NA,
              name_to_check = taxon_name
            ))
          }

          tryCatch({
            # Add small delay to respect API rate limits
            Sys.sleep(0.1)

            # Search for the taxon in WoRMS
            records <- wm_records_names(taxon_name, marine_only = FALSE)

            if (length(records) > 0 && !is.null(records[[1]]) && nrow(records[[1]]) > 0) {
              record <- records[[1]][1, ]

              return(data.frame(
                is_marine = ifelse(is.null(record$isMarine), NA, record$isMarine),
                is_brackish = ifelse(is.null(record$isBrackish), NA, record$isBrackish),
                is_freshwater = ifelse(is.null(record$isFreshwater), NA, record$isFreshwater),
                is_terrestrial = ifelse(is.null(record$isTerrestrial), NA, record$isTerrestrial),
                worms_rank = ifelse(is.null(record$rank), NA, record$rank),
                aphia_id = ifelse(is.null(record$AphiaID), NA, record$AphiaID),
                name_to_check = taxon_name
              ))

            } else {

              return(data.frame(
                is_marine = NA,
                is_brackish = NA,
                is_freshwater = NA,
                is_terrestrial = NA,
                worms_rank = NA,
                aphia_id = NA,
                name_to_check = taxon_name
              ))

            }
          },
          error = function(e) {
            message(paste("Error for:", taxon_name, "-", e$message))
            return(data.frame(
              is_marine = NA,
              is_brackish = NA,
              is_freshwater = NA,
              is_terrestrial = NA,
              worms_rank = NA,
              aphia_id = NA,
              name_to_check = taxon_name
            ))
          })
        }) |>
        bind_rows()


      gbif |>
        left_join(marinekey, by = "name_to_check") |>
        filter(is_marine == 1) |>
        select(-names(marinekey))



      }
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
        pred_within(st_as_text(simplecoastalmar$geoms)),
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
        filter(institutionCode == "iNaturalist")
      }
  ),

  tar_target(
    name = inatdataall,
    command = {
      gbifdataall |>
        filter(institutionCode == "iNaturalist")
      }
  ),

  tar_target(
    name = rinatdata,
    command = {
      bbox <- st_bbox(coastalmar)

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
          st_join(coastalmar, join = st_within)
      } else {
        inat_sf <- NULL
      }

      inat_sf
    }
  )
)


