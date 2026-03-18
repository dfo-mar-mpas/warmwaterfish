#R gbif functions

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
                      return(res$usageKey)
                             }



############################################################
# ---- 1. Get AquaMaps raster for species ----
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
# ---- 2. Metrics ----
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
# ---- 3. Iterative + checkpoint-safe loop ----
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


