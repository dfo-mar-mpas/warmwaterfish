#gbif credentials -----
Sys.setenv(
  GBIF_USER = "rystanley",
  GBIF_PWD  = "Homarus23",
  GBIF_EMAIL = "ryan.stanley@dfo-mpo.gc.ca"
)

#load libraries ---
pkgs <- c("tidyverse","sf","rnaturalearth",
          "terra","tidyterra","viridis",
          "targets","rgbif","purrr","qs",
          "rnaturalearth","aquamapsdata","DBI",
          "ggspatial","cli","rfishbase")

librarian::shelf(pkgs)

#load functions
source("code/00_functions.R")
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/trim_img_ws.R")

#projections ------
latlong <- 4326
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

############################################################
# ---- 0. Load data  ----
############################################################
dir("_targets/objects/")

  #shapefiles
  maritimes_sf <- tar_read(maritimes)
  coastal_sf <- tar_read(coastalmar)

  gbif_marine_data <- tar_read(gbifdataall)

  gbif_flat <- gbif_marine_data%>%
    st_drop_geometry

  #clean ws
  rm(gbif_marine_data);gc()

############################################################
# ---- 1. identify probable warm water species based on observation frequency ----
############################################################

  warm_sp_manual <- read.csv("data/GBIF_keys.csv")%>%
                    rbind(.,data.frame(Common.Name = c("Barracudas","Snakefishes"),
                                       Taxon.Name = c("Sphyraenidae","Trichinocephalus myops"),
                                       Taxon.Level = c("family","species"),
                                       GBIF.Key = c(NA,NA) # res <- name_backbone(name = "Trichinocephalus myops");res$usageKey
                                       ))%>%
                    rowwise()%>%
                    mutate(speciesKey = gbif_name(Taxon.Name,Taxon.Level))%>%
                    data.frame()

# -------------------
#identify rare species
# -------------------

  gbif_filtered <- gbif_flat %>%
    filter(phylum == "Chordata",
           !class %in% c("Myxini","Petromyzonti","Aves","Ascidiacea",
                         "Mammalia","Thaliacea","Testudines","Squamata"),#just keep sharks and fish. for some reason fish do not have a consistent class
           species != "Lycodichthys dearborni",
           occurrenceStatus == "PRESENT")

raresp1 <-  gbif_filtered %>%  #This is a misnomer from the pull antarctic species that was mis assigned to the region
          count(speciesKey, sort = TRUE)%>%
          filter(n<=30, #30 record threshold
                !is.na(speciesKey))%>%
          left_join(gbif_flat%>%
                      distinct(speciesKey,.keep_all=TRUE)%>%
                      dplyr::select(speciesKey,kingdom,phylum,class,order,family,genus,species))

#identify warm water species that were identified for the paper but removed based on the 'rarity' criterion above
missing_targets <- c(
                    setdiff(warm_sp_manual%>%filter(Taxon.Level=="genus")%>%pull(Taxon.Name),
                            raresp1%>%distinct(genus,.keep_all=TRUE)%>%pull(genus)),

                    setdiff(warm_sp_manual%>%filter(Taxon.Level=="family")%>%pull(Taxon.Name),
                            raresp1%>%distinct(family,.keep_all=TRUE)%>%pull(family)),

                    setdiff(warm_sp_manual%>%filter(Taxon.Level=="species")%>%pull(Taxon.Name),
                            raresp1%>%distinct(species,.keep_all=TRUE)%>%pull(species))
                  )%>%
                  setdiff(.,c("Aulostomus","Trichinocephalus myops"))#these are not in the data

#get the keys for the missing targets
    missing_target_speciesKey <- NULL
    for(i in missing_targets){

      temp <- warm_sp_manual%>%filter(Taxon.Name == i)

      temp2 <- gbif_flat%>%
               distinct(speciesKey,.keep_all=TRUE)

      if(temp$Taxon.Level == "genus"){temp3 <- temp2%>%
                                                filter(genus == temp$Taxon.Name)}

      if(temp$Taxon.Level == "family"){temp3 <- temp2%>%
                                        filter(family == temp$Taxon.Name)}

      if(temp$Taxon.Level == "species"){temp3 <- temp2%>%
                                        filter(species == temp$Taxon.Name)}

      missing_target_speciesKey <- c(missing_target_speciesKey,temp$speciesKey)

    }

  raresp <- c(raresp1%>%pull(speciesKey)%>%unique()%>%as.numeric,
              missing_target_speciesKey%>%as.numeric())

  #save the outputs
    write.csv(data.frame(speciesKey=raresp)%>%left_join(gbif_flat%>%
                                                   distinct(speciesKey,.keep_all=TRUE)%>%
                                                   dplyr::select(speciesKey,kingdom,phylum,class,order,family,genus,species)),
          file = "output/rare_species_30obs.csv",row.names=FALSE)


#### medium observation filter

    medsp1 <- gbif_filtered%>%  #This is a misnomer from the pull antarctic species that was mis assigned to the region
              count(speciesKey, sort = TRUE) %>%
              filter(n<=2000, #2000 record threshold
                     !is.na(speciesKey)) %>%
              pull(speciesKey) %>%
              as.numeric()%>%
              setdiff(.,raresp)

    med_species <- medsp1%>%
                   data.frame(speciesKey=.)%>%
                   left_join(gbif_filtered%>%
                                distinct(speciesKey,.keep_all=TRUE)%>%
                                dplyr::select(speciesKey,kingdom,phylum,class,order,family,genus,species))%>%
                   left_join(gbif_filtered%>%  #This is a misnomer from the pull antarctic species that was mis assigned to the region
                               count(speciesKey, sort = TRUE)%>%
                               dplyr::select(speciesKey,n))

#map out the study range------

    basemap <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                       dplyr::select(name_en,geometry)%>%
                       st_as_sf()%>%
                       st_union()%>%
                       st_transform(latlong)%>%
                       st_as_sf()%>%
                       mutate(country="Canada"),
                     ne_states(country = "United States of America",returnclass = "sf")%>%
                       dplyr::select(name_en,geometry)%>%
                       st_as_sf()%>%
                       st_union()%>%
                       st_transform(latlong)%>%
                       st_as_sf()%>%
                       mutate(country="USA"))

    #for rgbif note a range to acruire data. In this case the Northwest Atlantic (Caribean to Arctic) to assess 'climate migrants'
    filter_range <- st_bbox(c(xmin = -100, xmax = -40, ymin = 13, ymax = 85), crs = 4326)

    #create a projection for a spherical map
    center_pt <- filter_range%>%
      st_as_sfc()%>%
      st_transform(latlong)%>%
      st_centroid()

    lon0 <- st_coordinates(center_pt)[1]
    lat0 <- st_coordinates(center_pt)[2]

    globe_crs <- sprintf("+proj=ortho +lat_0=%s +lon_0=%s",lat0, lon0)

    #download the world globe basemap
    world_globe <- ne_countries(
      scale = "medium",
      returnclass = "sf"
    ) %>%
      st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>%
      st_transform(globe_crs)

    globe <- ne_countries(
      scale = "medium",
      returnclass = "sf"
    )%>%st_transform(latlong)

    #define the box denoting the study region you want to highlight
    global_box <- filter_range%>%
      st_as_sfc() %>%
      st_transform(globe_crs)

    #make a circle to wrap the globe plot
    globe_circle <- st_sfc(
      st_buffer(
        st_point(c(0, 0)),   # center of orthographic projection
        dist =  6378137  # meters
      ),
      crs = globe_crs
    )

    #crudgy way to make it so that the oceans are white in the plot
    globe_disc <- st_sfc(
      st_point(c(0, 0)),  # center in projected coords
      crs = globe_crs
    ) %>%
      st_buffer(dist = 6378137) %>%   # Earth radius in meters
      st_as_sf()

    #MAKE INSET PLOT
    global_inset <- ggplot()+
                geom_sf(data = globe_disc, fill = "white", colour = "black", linewidth = 0.4)+
                geom_sf(data = world_globe,colour = "grey20",linewidth = 0.2)+
                geom_sf(data = world_globe%>%filter(formal_en == "Canada"), fill = "grey60",colour = "grey20",linewidth = 0.2)+
                geom_sf(data = global_box,fill = NA,colour = "black",linewidth = 0.3)+
                geom_sf(data = globe_circle,fill = NA,colour = "grey30",linewidth = 0.4)+
                coord_sf(crs = globe_crs) +
                theme_void() +
                theme(panel.background = element_rect(fill = NA, colour = NA),
                  plot.background  = element_rect(fill = NA, colour = NA))

    #make primary plots
    atlantic_range <- coastal_sf%>%
                      st_bbox()%>%
                      st_as_sfc()%>%
                      st_transform(utm)%>%
                      st_buffer(20)%>%
                      st_transform(latlong)%>%
                      st_bbox()

      p1 <- ggplot()+
        geom_sf(data=globe%>%st_transform(latlong))+
        geom_sf(data = world_globe%>%st_transform(latlong)%>%filter(formal_en == "Canada"),
          fill = "grey60")+
        geom_sf(data=atlantic_range%>%st_as_sfc(),fill=NA)+
        coord_sf(xlim=filter_range[c(1,3)],ylim=filter_range[c(2,4)],expand=0)+
        theme_bw()

      p1_inset <- p1 +
        inset_element(
          global_inset,
          left = 0.65,   # adjust position
          bottom = 0.05,
          right = 0.98,
          top = 0.4,
          align_to = "panel"
        )

      p2 <- ggplot()+
        geom_sf(data=coastal_sf%>%st_transform(latlong),fill="cornflowerblue")+
        geom_sf(data=basemap)+
        geom_sf(data=basemap%>%filter(country=="Canada"),fill = "grey60")+
        coord_sf(xlim=atlantic_range[c(1,3)],ylim=atlantic_range[c(2,4)],expand=0)+
        theme_bw()+
        annotation_scale(location = "br")

      #make a combination plot
      combo <- p1_inset+p2+plot_layout(ncol=2)

      ggsave("output/mapping_region.png",combo,height=6,width=8,units="in",dpi=300)
      trim_img_ws("output/mapping_region.png") #helps with the need to finagle the dimensions above

############################################################
# ---- 1. Download GBIF data for rare species ----
############################################################

      #set bounding box
      bbox <- filter_range |>
        st_as_sfc() |>
        st_as_text()

      #download all of it together
      # #establish download key with gbif
      download_key <- occ_download(
        pred_within(bbox), #for the box mapped above
        pred("hasCoordinate", TRUE),
        pred("hasGeospatialIssue", FALSE),
        pred_in("speciesKey", c(raresp,medsp1)), #this is only for the records specified above
        format = "SIMPLE_CSV",
        user = Sys.getenv("GBIF_USER"),
        pwd = Sys.getenv("GBIF_PWD"),
        email = Sys.getenv("GBIF_EMAIL")
      )

      # check status of download key
      occ_download_meta(download_key)

      # when status == "SUCCEEDED"
      dat <- occ_download_get(download_key, path = "data/gbif_download", overwrite = FALSE)
      save(dat,file="data/gbif_dat_get_all.RData")

#load the download
# import as data frame
load("data/gbif_dat_get_all.RData")
gbif_df <- occ_download_import(dat)

#missing one species
sennet <- (rgbif::occ_data(
  scientificName = "Sphyraena borealis",
  hasCoordinate = TRUE,
  limit = 50000,
  decimalLongitude = "-100,-10"  # Western Atlantic filter only
    )$data)%>%
  mutate(class="Teleostei")

sennet <- sennet%>%
         dplyr::select(intersect(names(gbif_df),names(sennet)))

gbif_df <- gbif_df%>%
           dplyr::select(intersect(names(gbif_df),names(sennet)))%>%
           rbind(.,sennet)


#compute summary stats for each species
species_stats <- gbif_df %>%
  filter(phylum == "Chordata",
         !class %in% c("Myxini","Petromyzonti","Aves","Ascidiacea",
                       "Mammalia","Thaliacea","Testudines","Squamata"),#just keep sharks and fish. for some reason fish do not have a consistent class
         species != "Lycodichthys dearborni",
         occurrenceStatus == "PRESENT")%>%
  group_by(species) %>%
  summarise(
    n = n(),
    q95 = quantile(decimalLatitude, 0.95, na.rm=TRUE),
    q99 = quantile(decimalLatitude, 0.99, na.rm=TRUE),
    median_lat = median(decimalLatitude,na.rm=TRUE),
    mean_lat = median(decimalLongitude,na.rm=T)
  )%>%
  mutate(
    warm_flag = case_when(
      q95 < 38 ~ "tropical",
      q95 < 41 ~ "southern_range_edge",
      q95 < 43.2 ~ "possible_warm_water",
      TRUE ~ "northern_species"
    )
  )%>%
  left_join(.,gbif_df%>%
              distinct(species,.keep_all=TRUE)%>%
              dplyr::select(kingdom,phylum,class,order,genus,family,species,speciesKey))%>%
  left_join(.,data.frame(rbind(data.frame(speciesKey=raresp,selection="rare"),
                               data.frame(speciesKey=medsp1,selection="less than 2000 obs"))))

#save the output --
write.csv(species_stats,"output/species_stats.csv",row.names=FALSE)

############################################################
# ----  Aquamaps integration  ----
############################################################

# ----  One-time setup (run once) ----
# Install AquaMaps database (~2GB download, ~10GB unpacked)
#download_db()

# Set database type
default_db("sqlite")

# Connect to AquaMaps database
con <- con_am("sqlite")

############################################################
# ---- flag anythign missing from the table ----
############################################################

#some species are missing from the GBIF pull from the warm water species table in the figure
missing_targets <- c(setdiff(warm_sp_manual%>%filter(Taxon.Level=="genus")%>%pull(Taxon.Name),
                            species_stats%>%distinct(genus,.keep_all=TRUE)%>%pull(genus)),

                    setdiff(warm_sp_manual%>%filter(Taxon.Level=="family")%>%pull(Taxon.Name),
                            species_stats%>%distinct(family,.keep_all=TRUE)%>%pull(family)),

                    setdiff(warm_sp_manual%>%filter(Taxon.Level=="species")%>%pull(Taxon.Name),
                            species_stats%>%distinct(species,.keep_all=TRUE)%>%pull(species))
                    )%>%setdiff(.,c("Aulostomus","Trichinocephalus myops")) #nothing from the tropical paper table is missing


############################################################
# ---- Initial aquamaps extractions  ----
############################################################

#Run the extractions from Aquamaps. note any NA summaries means that the distribution from aquamps for that species is not within the bounding box of this analysis (western atlantic, north of the equator)
  aquamaps_process <- run_aquamaps_analysis(
                      species_vec = species_stats$species,
                      coastal_sf = coastal_sf,
                      prob_cutoff = threshold,
                      save_file="output/aquamaps_results.csv" #this checks to see if it has already been done
                    )

  aquamaps_results <- read.csv("output/aquamaps_results.csv")%>%
                      filter(species%in%species_stats$species)

############################################################
# ---- key out problem species and fix them  ----
############################################################

problem_species <- aquamaps_results%>%filter(is.na(weighted_lat)|is.na(q95_lat)|is.na(coastal_overlap))%>%pull(species)

  #run checks
    #1) does the name work with aquamaps, is there an accepted name to replace the synonym - aquamaps results can be reassessed using the accepted name in this case
    #2) if name is good, does aquamaps have a map for it
    #3) if yes to 1) and 2) where is the species based on the fishbase FAO region 21/31 being the study area and if no FAO is the 'country' recorded for that species

species_problem_table <- resolve_species_vector(problem_species)

############################################################
# ---- Knit it all together  ----
############################################################


aquamaps_results <- read.csv("output/aquamaps_results.csv")%>%
                    filter(species %in% c(species_stats$species,missing_targets%>%pull(species2)))%>%
                    mutate(warm_flag_am = case_when(
                      q95_lat < 38 ~ "tropical",
                      q95_lat < 41 ~ "southern_range_edge",
                      q95_lat < 43.2 ~ "possible_warm_water",
                      is.na(q95_lat) ~ NA,
                      TRUE ~ "northern_species"
                    ))

warm_water_summary <- species_stats%>%
                      left_join(.,aquamaps_results)%>%
                      mutate()


aquamaps_results_rare <- read.csv("output/aquamaps_results.csv")%>% #reload and trim
                         filter(species %in% tropical_sp)%>%
                         mutate(warm_flag_am = case_when(
                           q95_lat < 38 ~ "tropical",
                           q95_lat < 41 ~ "southern_range_edge",
                           q95_lat < 43.2 ~ "possible_warm_water",
                           is.na(q95_lat) ~ NA,
                           TRUE ~ "northern_species"
                         ))%>%mutate(species_selection="rare")

med_tropical_species <- species_stats%>%filter(species_selection=="medium",warm_flag == "tropical")

aquamaps_results_med <- read.csv("output/aquamaps_results.csv")%>% #reload and trim
                          filter(species %in% med_species$species)%>%
                          mutate(warm_flag_am = case_when(
                            q95_lat < 38 ~ "tropical",
                            q95_lat < 41 ~ "southern_range_edge",
                            q95_lat < 43.2 ~ "possible_warm_water",
                            is.na(q95_lat) ~ NA,
                            TRUE ~ "northern_species"
                          ))%>%mutate(species_selection="medium")




ggplot(warm_water_summary,aes(q95,q95_lat))+
  geom_point()+
  theme_bw()+
  geom_hline(yintercept=43,linetype=2)+
  stat_smooth(method="lm")+
  labs(x="GBIF 95th percentile °lat",
       y="Aquamaps weighted mean °lat")

###########################################################
# ---- raster plots tropical species matching to Hunter's montage (Figure 1) ----
############################################################
Figure1_df <- data.frame(key=LETTERS[1:10],
                          species=c("Balistes vetula","Pristigenys alta","Acanthostracion quadricornis",
                                    "Hyporthodus niveatus","Dactylopterus volitans","Heteropriacanthus cruentatus",
                                    "Chaetodon ocellatus","Mullus auratus","Seriola zonata","Hippocampus erectus"))

# -----------------------------
# Get aquamap rasters
# -----------------------------
species_rasters <- Figure1_df %>%
  mutate(rast = map(species, get_aquamaps_raster)) %>%
  filter(!map_lgl(rast, is.null))

# -----------------------------
# Build common raster template
# -----------------------------
# AquaMaps is ~0.5° grid → enforce that
# template <- rast(
#   xmin = filter_range$xmin,
#   xmax = filter_range$xmax,
#   ymin = filter_range$ymin,
#   ymax = filter_range$ymax,
#   resolution = 0.5,
#   crs = "EPSG:4326"
# )

template <- rast(
  xmin = -100,
  xmax = -30,
  ymin = -60,
  ymax = 85,
  resolution = 0.5,
  crs = "EPSG:4326"
)

# -----------------------------
# Align all rasters to template
# -----------------------------
species_rasters <- species_rasters %>%
  mutate(
    rast = map(rast, ~{
      r <- .x

      # project if needed
      if (!raster::compareCRS(r, template)) {
        r <- project(r, template)
      }

      # snap to common grid
      r <- resample(r, template, method = "bilinear")

      r
    })
  )

# -----------------------------
# Stack, crop, apply threshold
# -----------------------------

threshold <- 0.5

rast_stack <- rast(species_rasters$rast)
names(rast_stack) <- species_rasters$key
#rast_stack <- rast_stack%>%crop(.,filter_range%>%st_as_sfc()%>%vect())
rast_stack[rast_stack < threshold] <- NA

# -----------------------------
# Weighted latitude
# -----------------------------
weighted_lats <- map_dfr(1:nlyr(rast_stack), function(i){

  df <- as.data.frame(rast_stack[[i]], xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "value"

  tibble(
    lyr = names(rast_stack)[i],
    pw_lat = weighted.mean(df$y, df$value)
  )
})


# -----------------------------
# Plot the rasters as a facet
# -----------------------------
raster_facet <- ggplot() +
                geom_spatraster(data = rast_stack) +
                geom_sf(data = globe) +
                geom_sf(data = globe%>%filter(formal_en == "Canada"), fill = "grey60") +
                geom_hline(
                  data = weighted_lats,
                  aes(yintercept = pw_lat),
                  linetype = "dashed",
                  colour = "red",
                  linewidth = 0.5
                ) +
                geom_sf(data=atlantic_range%>%st_as_sfc(),fill=NA)+
                scale_fill_viridis_c(option = "viridis",limits = c(0, 1),na.value = NA) +
                facet_wrap(~lyr, ncol = 5) +
                coord_sf(xlim = c(filter_range[1],-55),ylim = c(filter_range[2],50),expand = 0) +
                scale_x_continuous(breaks = seq(-90, -50, by = 10))+
                theme_bw() +
                theme(
                  strip.background = element_rect(fill = "white"),
                  strip.text = element_text(size = 10),
                  axis.text.y = element_text(size=6),
                  axis.text.x = element_text(angle = 30, hjust = 1,size=6),
                  panel.grid = element_blank()
                ) +
                labs(fill = "Probability")

ggsave("output/raster_facet.png",raster_facet,height=6,width=8.5,units="in",dpi=300)
trim_img_ws("output/raster_facet.png")
