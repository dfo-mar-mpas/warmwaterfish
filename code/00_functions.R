############################################################
# ---- GBIF data processing functions ----
############################################################

get_lat_stats <- function(species,
                          lim = 50000,
                          lat_thresh = 43.2,
                          cache_dir = "data/occ_cache"){

  # sanitize filename
  file_name <- paste0(gsub(" ", "_", species), ".qs")
  file_path <- file.path(cache_dir, file_name)

  # load cached data if it exists
  if(file.exists(file_path)){

    occ <- qs::qread(file_path)

  } else {

    occ <- tryCatch(
        rgbif::occ_data(
        scientificName = species,
        hasCoordinate = TRUE,
        limit = lim,
        decimalLongitude = "-100,-10"  # Western Atlantic filter only
      )$data,
      error = function(e) return(NULL)
    )

    if(is.null(occ) || nrow(occ) == 0) return(NULL)

    # create cache directory if needed
    dir.create(cache_dir, showWarnings = FALSE)

    # save raw occurrence data
    qs::qsave(occ, file_path)
  }

  # quick data cleaning
  occ <- occ %>%
    dplyr::filter(
      !(decimalLongitude == 0 & decimalLatitude == 0),
      decimalLongitude >= -100,
      decimalLongitude <= -30
    )

  if(nrow(occ) == 0) return(NULL)

  lat <- occ$decimalLatitude
  lat <- sample(lat, min(length(lat), 20000)) #random sub sample to speed up the processing

  data.frame(
    species = species,
    n_records = length(lat),
    mean_lat = mean(lat, na.rm = TRUE),
    median_lat = median(lat, na.rm = TRUE),
    max_lat = max(lat, na.rm = TRUE),
    min_lat = min(lat, na.rm = TRUE),
    quant_90 = as.numeric(quantile(lat, 0.9, na.rm = TRUE)),
    prop_above_thresh = mean(lat >= lat_thresh, na.rm = TRUE)
  )
}

get_lat_stats_triplicate <- function(species,
                                     lim = 200,
                                     lat_thresh = 43.2,
                                     reps = 3){

  run_query <- function(){

    occ <- tryCatch(
      rgbif::occ_data(
        scientificName = species,
        hasCoordinate = TRUE,
        limit = lim,
        decimalLongitude = "-100,-10"
      )$data,
      error = function(e) return(NULL)
    )

    if(is.null(occ) || nrow(occ) == 0) return(NULL)

    occ <- occ %>%
      filter(
        !(decimalLongitude == 0 & decimalLatitude == 0),
        decimalLongitude >= -100,
        decimalLongitude <= -30
      )

    if(nrow(occ) == 0) return(NULL)

    lat <- occ$decimalLatitude

    tibble(
      q90 = quantile(lat,0.9,na.rm=TRUE),
      q95 = quantile(lat,0.95,na.rm=TRUE),
      q99 = quantile(lat,0.99,na.rm=TRUE),
      prop_above = mean(lat >= lat_thresh,na.rm=TRUE)
    )
  }

  reps_out <- map_dfr(1:reps, ~run_query())

  if(nrow(reps_out) == 0) return(NULL)

  tibble(
    species = species,
    q90_mean = mean(reps_out$q90),
    q90_sd   = sd(reps_out$q90),
    q95_mean = mean(reps_out$q95),
    q95_sd   = sd(reps_out$q95),
    q99_mean = mean(reps_out$q99),
    q99_sd   = sd(reps_out$q99),
    prop_above_mean = mean(reps_out$prop_above),
    prop_above_sd   = sd(reps_out$prop_above)
  )
}

gbif_name <- function(Taxon.Name,Taxon.Level){
                      if(Taxon.Level == "species"){res <- name_backbone(name = Taxon.Name)}
                      if(Taxon.Level == "genus"){res <- name_backbone(genus = Taxon.Name)}
                      if(Taxon.Level == "family"){res <- name_backbone(family = Taxon.Name)}
                      return(res$usageKey)}

############################################################
# ---- Get AquaMaps raster for species ----
############################################################

get_aquamaps_cells <- function(species, bound_box = NULL){

  # ---- Default Western Atlantic bounding box ----
  if(is.null(bound_box)){
    bound_box <- st_bbox(c(
      xmin = -100,
      xmax = -40,
      ymin = 13,
      ymax = 85
    ), crs = st_crs(4326))
  }

  # ---- Get species ID (fuzzy match) ----
  sp_match <- tryCatch(
    am_search_fuzzy(species),
    error = function(e) return(NULL)
  )

  if(is.null(sp_match) || nrow(sp_match) == 0) return(NULL)
  sp_id <- sp_match$key[1]

  # ---- Pull raster ----
  ras <- tryCatch(
    am_raster(sp_id),
    error = function(e) return(NULL)
  )
  if(is.null(ras)) return(NULL)

  # ---- Convert to terra raster and project ----
  r <- terra::rast(ras)
  r <- terra::project(r, "EPSG:4326")

  # ---- Crop to bounding box ----
  bb_vect <- terra::vect(st_as_sfc(bound_box))
  r_crop <- tryCatch(
    terra::crop(r, bb_vect),
    error = function(e) return(NULL)
  )

  # ---- Check for empty crop ----
  if(is.null(r_crop)) return(NULL)
  if(terra::ncell(r_crop) == 0) return(NULL)  # no overlapping cells

  # ---- Convert to data.frame ----
  df <- terra::as.data.frame(r_crop, xy = TRUE, na.rm = TRUE)

  if(nrow(df) == 0) return(NULL)  # safety check
  names(df) <- c("lon", "lat", "probability")

  return(df)
}

get_aquamaps_raster <- function(species){

  sp_match <- tryCatch(
    am_search_fuzzy(species),
    error = function(e) return(NULL)
  )

  if(is.null(sp_match) || nrow(sp_match) == 0) return(NULL)

  sp_id <- sp_match$key[1]

  ras <- tryCatch(
    am_raster(sp_id),
    error = function(e) return(NULL)
  )

  if(is.null(ras)) return(NULL)

  ras <- terra::rast(ras)

  return(ras)
}

############################################################
# ---- Metrics ----
############################################################

calc_weighted_lat <- function(df, prob_cutoff = 0){

  df <- df %>% filter(probability > prob_cutoff)
  if(nrow(df) == 0) return(NA_real_)

  sum(df$lat * df$probability, na.rm = TRUE) /
    sum(df$probability, na.rm = TRUE)
}

calc_weighted_q95 <- function(df, prob_cutoff = 0){

  df <- df %>%
    filter(probability > prob_cutoff) %>%
    arrange(lat)

  if(nrow(df) == 0) return(NA_real_)

  df %>%
    mutate(cumprob = cumsum(probability) / sum(probability)) %>%
    filter(cumprob >= 0.95) %>%
    slice(1) %>%
    pull(lat)
}

calc_coastal_overlap <- function(df, coastal_sf, prob_cutoff = 0){

  df <- df %>% filter(probability > prob_cutoff)
  if(nrow(df) == 0) return(NA_real_)

  pts <- st_as_sf(df, coords = c("lon","lat"), crs = 4326)

  inside <- st_intersects(pts, coastal_sf, sparse = FALSE)

  sum(df$probability[inside], na.rm = TRUE) /
    sum(df$probability, na.rm = TRUE)
}

############################################################
# ---- clean up taxonomy to link to aqumaps ----
############################################################

resolve_aquamaps_species <- function(species, bound_box = NULL, atlantic_vect = NULL){
  library(terra)
  library(sf)
  library(rfishbase)
  library(rnaturalearth)

  # -----------------------------
  # Default Western Atlantic bounding box
  # -----------------------------
  if(is.null(bound_box)){
    bound_box <- st_bbox(c(
      xmin = -100, xmax = -40,
      ymin = 13,  ymax = 85
    ), crs = st_crs(4326))
  } #note the southwestern corner of this box does include part of the pacific


  # -----------------------------
  # FAO reference table
  # -----------------------------
  fao_ref <- data.frame(
    AreaCode = c(18,21,27,31,34,37,41,47,48,51,57,58,61,67,71,77,81,87,88),
    AreaName = c(
      "Arctic Sea","Atlantic, Northwest","Atlantic, Northeast","Atlantic, Western Central",
      "Atlantic, Eastern Central","Mediterranean and Black Sea","Atlantic, Southwest",
      "Atlantic, Southeast","Atlantic, Antarctic","Indian Ocean, Western","Indian Ocean, Eastern",
      "Indian Ocean, Antarctic and Southern","Pacific, Northwest","Pacific, Northeast",
      "Pacific, Western Central","Pacific, Eastern Central","Pacific, Southwest",
      "Pacific, Southeast","Pacific, Antarctic"
    )
  )
  study_region_codes <- c(21,31) # Western Atlantic codes

  # -----------------------------
  # Helper: try AquaMaps
  # -----------------------------
  try_aquamaps <- function(name){
    res <- tryCatch(am_search_fuzzy(name), error = function(e) NULL)
    if(!is.null(res) && nrow(res) > 0){
      terms_str <- res$terms[1]
      sci_name <- NA_character_
      if(!is.na(terms_str)){
        split_terms <- strsplit(terms_str, " ")[[1]]
        if(length(split_terms) >= 2) sci_name <- paste(split_terms[1:2], collapse = " ")
      }
      return(list(name = sci_name, key = res$key[1]))
    }
    return(NULL)
  }

  #do the FAO check - fishbase 3rd way to match to the region
  fao_info <- tryCatch(faoareas(species), error = function(e) NULL)%>%suppressMessages()

  if(!is.null(fao_info) && nrow(fao_info) > 0 && "AreaCode" %in% colnames(fao_info)){

    fao_region <- paste(fao_ref$AreaName[fao_ref$AreaCode %in% fao_info$AreaCode], collapse = "_")
    in_study_region <- any(fao_info$AreaCode %in% study_region_codes)

  } else {
    fao_region <- NA_character_
    in_study_region = NA

  }

  # -----------------------------
  # 1. Try original name in AquaMaps
  # -----------------------------
  out <- try_aquamaps(species)
  fao_region <- NA_character_
  in_study_region <- NA

  if(!is.null(out)){
    ras <- tryCatch(am_raster(out$key), error = function(e) NULL)

    if(!is.null(ras)){
      ras <- terra::rast(ras) %>% terra::project("EPSG:4326")
      bb_vect <- bound_box %>% st_as_sfc() %>% terra::vect()
      ras_crop <- tryCatch(terra::crop(ras, bb_vect), error = function(e) NULL)

      if (!is.null(ras_crop)) {
        ras_mask <- tryCatch(terra::mask(ras_crop, atlantic_vect), error = function(e) NULL)
      } else {
        ras_mask <- NULL
      }

      if(!is.null(ras_crop) && !is.na(sum(terra::values(ras_mask), na.rm = TRUE)) &&
         sum(terra::values(ras_mask), na.rm = TRUE) > 0){
        # Raster intersects bbox → success
        return(data.frame(
          input_name      = species,
          matched_name    = out$name,
          speciesKey      = out$key,
          status          = "aquamaps_match",
          FAO_region      = fao_region,
          in_study_region = in_study_region,
          country         = NA_character_
        ))
      } else {
          return(data.frame(
          input_name      = species,
          matched_name    = out$name,
          speciesKey      = out$key,
          status          = "aquamaps_outside_bbox",
          FAO_region      = fao_region,
          in_study_region = in_study_region,
          country         = NA_character_
        ))
      }
    } else {
      # No raster → check FAO

      return(data.frame(
        input_name      = species,
        matched_name    = out$name,
        speciesKey      = out$key,
        status          = "aquamaps_no_raster",
        FAO_region      = fao_region,
        in_study_region = in_study_region,
        country         = NA_character_
      ))
    }
  }

  # -----------------------------
  # 2. FishBase synonyms
  # -----------------------------
  fb_syn <- tryCatch(synonyms(species), error = function(e) NULL)
  candidate_names <- species
  if(!is.null(fb_syn) && nrow(fb_syn) > 0){
    candidate_names <- unique(c(fb_syn$Species, fb_syn$Synonym))
    candidate_names <- candidate_names[!is.na(candidate_names)]
  }

  # -----------------------------
  # 3. Try all FishBase names in AquaMaps
  # -----------------------------
  for(nm in candidate_names){
    out <- try_aquamaps(nm)
    if(!is.null(out)){
      return(data.frame(
        input_name      = species,
        matched_name    = out$name,
        speciesKey      = out$key,
        status          = "fishbase_synonym_corrected",
        FAO_region      = fao_region,
        in_study_region = in_study_region,
        country         = NA_character_
      ))
    }
  }

  # -----------------------------
  # 4. GBIF fallback
  # -----------------------------
  gbif_match <- tryCatch(name_backbone(name = species), error = function(e) NULL)
  if(!is.null(gbif_match) && !is.na(gbif_match$scientificName)){
    out <- try_aquamaps(sub("^((\\S+\\s+\\S+)).*", "\\1",gbif_match$scientificName)) #remove the taxonomic authority that gbif returns but is not used by aquamaps
    if(!is.null(out)){
      return(data.frame(
        input_name      = species,
        matched_name    = out$name,
        speciesKey      = out$key,
        status          = "gbif_corrected",
        FAO_region      = fao_region,
        in_study_region = in_study_region,
        country         = NA_character_
      ))
    }
  }

  # -----------------------------
  # 5. FishBase FAO / country fallback
  # -----------------------------
  fao_info <- tryCatch(faoareas(species), error = function(e) NULL)
  fao_region <- NA_character_
  in_study_region <- NA
  country_vec <- NA_character_

  if(!is.null(fao_info) && nrow(fao_info) > 0 && "AreaCode" %in% colnames(fao_info)){
    fao_region <- paste(fao_ref$AreaName[fao_ref$AreaCode %in% fao_info$AreaCode], collapse = "_")
    in_study_region <- any(fao_info$AreaCode %in% study_region_codes)
  } else {
    fb_ctry <- tryCatch(country(species), error = function(e) NULL)
    if(!is.null(fb_ctry) && nrow(fb_ctry) > 0 && "Country" %in% colnames(fb_ctry)){
      countries <- unique(fb_ctry$Country)
      countries <- countries[!is.na(countries)]
      if(length(countries) > 0) country_vec <- paste(countries, collapse = "_")
    }
  }

  status2 <- if (is.na(fb_syn$Species[1])) {
    "fishbase_not_assessed"
  } else if (species == fb_syn$Species[1]) { #catch scenarios where the synonym is the same as species - this would be in aquamaps shows 'accepted' but no map produced
    "aquamaps_no_raster"
  } else {
    "fishbase_not_assessed"
  }

  return(data.frame(
    input_name      = species,
    matched_name    = species,
    speciesKey      = NA,
    status          = status2,
    FAO_region      = fao_region,
    in_study_region = in_study_region,
    country         = country_vec
  ))
}

resolve_species_vector <- function(species_vec, bound_box = NULL, atlantic_vect = NULL){

  n <- length(species_vec)
  pb <- cli_progress_bar("Resolving species", total = n)

  all_results <- list()

  for(i in seq_along(species_vec)){
    sp <- species_vec[i]

    # Show which species is being processed
    message("Processing: ", sp)

    res <- resolve_aquamaps_species(sp, bound_box = bound_box, atlantic_vect=atlantic_vect)
    all_results[[i]] <- res

    # Update progress bar
    cli_progress_update(pb, set = i)
  }

  cli_progress_done(pb)

  # Combine all results into a single data frame
  species_table <- bind_rows(all_results)
  return(species_table)
}

fishbase_syn_check <- function(species) {
  fb_syn <- tryCatch(synonyms(species), error = function(e) NULL)

  if (is.null(fb_syn) || nrow(fb_syn) == 0 || is.na(fb_syn$Species[1])) {
    return(NA_character_)
  }

  fb_syn$Species[1]
} #just a double check because this doesn't quite work right in resolve aquamaps species.

aquamaps_check <- function(name){
  res <- tryCatch(am_search_fuzzy(name), error = function(e) NULL)
  if(!is.null(res) && nrow(res) > 0){
    terms_str <- res$terms[1]
    sci_name <- NA_character_
    if(!is.na(terms_str)){
      split_terms <- strsplit(terms_str, " ")[[1]]
      if(length(split_terms) >= 2) sci_name <- paste(split_terms[1:2], collapse = " ")
    }
    return(data.frame(name = sci_name, key = res$key[1]))
  }
  return(data.frame(name=NA,key=NA))
}

############################################################
# ---- Iterative + checkpoint-safe loop ----
############################################################

run_aquamaps_analysis <- function(species_vec,
                                  coastal_sf,
                                  prob_cutoff = 0.2,
                                  save_file = "aquamaps_results.csv"){

  results <- list()

  # Resume if file exists
  if(file.exists(save_file)){
    existing <- read_csv(save_file, show_col_types = FALSE)
    done_species <- existing$species
    species_vec <- setdiff(species_vec, done_species)
    results <- split(existing, seq(nrow(existing)))
  }

  for(i in seq_along(species_vec)){

    sp <- species_vec[i]
    message("Processing ", i, "/", length(species_vec), ": ", sp)

    df <- get_aquamaps_cells(sp)

    if(is.null(df)){
      res <- tibble(
        species = sp,
        weighted_lat = NA,
        q95_lat = NA,
        coastal_overlap = NA,
        n_cells = 0
      )
    } else {
      res <- tibble(
        species = sp,
        weighted_lat = calc_weighted_lat(df, prob_cutoff),
        q95_lat = calc_weighted_q95(df, prob_cutoff),
        coastal_overlap = calc_coastal_overlap(df, coastal_sf, prob_cutoff),
        n_cells = nrow(df)
      )
    }

    results[[length(results) + 1]] <- res

    # ---- Save progress every iteration ----
    write_csv(bind_rows(results), save_file)
  }

  bind_rows(results)
}


############################################################
# ---- Build an Atlantic mask ----
############################################################

# oceans <- ne_download(scale = 110, type = "geography_marine_polys", category = "physical", returnclass = "sf")
# land <- ne_countries(scale = 110, returnclass = "sf") %>% st_transform(4326)
#
# pacific <- oceans %>%
#   dplyr::filter(grepl("Pacific", name) )%>%
#   st_transform(4326)%>%
#   st_make_valid()%>%
#   st_union()
#
# land_mask <- land %>%
#   filter(admin %in% c(
#     "Mexico", "Belize", "Guatemala", "Honduras",
#     "El Salvador", "Nicaragua", "Costa Rica", "Panama"
#   )) %>%
#   st_union() %>%
#   st_make_valid()
#
# bb_sf <- st_as_sfc(bound_box)
#
# atlantic_bb <- bb_sf %>%
#   st_difference(pacific) %>%
#   st_difference(land_mask) %>%
#   st_make_valid()
#
# # Only keep largest contiguous part
# atlantic_bb <- st_cast(atlantic_bb, "POLYGON")
# atlantic_bb <- atlantic_bb[which.max(st_area(atlantic_bb)), ]
# atlantic_bb <- st_buffer(atlantic_bb, 0)   # fix tiny topology issues
#
# #this process does create some boundary effects (particulary in the south)
# south_fix <- bound_box
# south_fix[1] <- -85
# south_fix <- south_fix%>%st_as_sfc()
#
# atlantic_bb <- atlantic_bb%>%
#                st_union(.,south_fix)%>%
#                st_make_valid()
#
# write_sf(atlantic_bb,"data/Atlantic_bounding.shp")


