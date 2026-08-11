### Host preference analysis
### 2026.08.11
### Oliver G Spacey

# pre-amble ---------------------------------------------------------------
# clear environment
rm(list = ls())

# load packages
library(tidyverse) # data wrangling
library(ggplot2)   # data visualisation
library(terra)     # raster data
library(rgbif)      # GBIF occurrence data

# import data -------------------------------------------------------------
# All input files are stored in the project directory alongside this script.
# load host preference data
host_pref_df <- read_csv(
  "host_preference_data.csv",
  show_col_types = FALSE
)

# load list of hosts from Barney, Hawskworth and Geils, 1998
host_list_df <- read_csv(
  "viscum_album_subsp_album_hosts_current_names.csv",
  show_col_types = FALSE
)

# load koppen_geiger raster
kg_raster <- terra::rast("koppen_geiger_0p1.tif")

# load and parse the Köppen-Geiger legend
koppen_geiger_legend <- readLines("koppen_geiger_legend.txt")

koppen_geiger_legend_df <- tibble(
  legend_line = koppen_geiger_legend
) %>%
  filter(str_detect(legend_line, "^\\s*[0-9]+:")) %>%
  transmute(
    value = as.integer(str_extract(legend_line, "^[[:space:]]*[0-9]+")),
    koppen_geiger = str_match(legend_line, ":\\s+([A-Za-z]{2,3})\\s")[, 2],
    rgb = str_match(legend_line, "\\[([^]]+)\\]")[, 2]
  ) %>%
  mutate(rgb = str_squish(rgb))

# load trait data for European trees
eur_trait_df <- read_csv(
  "European_IW_data.csv",
  show_col_types = FALSE
)

# load trait data for North American trees
nam_trait_df <- read_csv(
  "N_America_IW_data.csv",
  show_col_types = FALSE
)

# wrangle data -------------------------------------------------------------
# Combine the regional wood-trait datasets and retain their source.
trait_dfs <- list(
  European = eur_trait_df,
  `North American` = nam_trait_df
)

# Some binary trait columns are imported with different types between files
# (for example, character in one file and numeric/logical in the other).
# Harmonise only those conflicting columns before binding the rows.
trait_columns <- names(eur_trait_df)
conflicting_trait_columns <- trait_columns[
  vapply(
    trait_columns,
    function(column) {
      length(unique(vapply(trait_dfs, function(df) typeof(df[[column]]), character(1)))) > 1
    },
    logical(1)
  )
]

trait_dfs <- lapply(
  trait_dfs,
  function(df) mutate(df, across(all_of(conflicting_trait_columns), as.character))
)

# combine datasets
wood_trait_df <- bind_rows(trait_dfs, .id = "dataset")

# Find obscure strings and remove
# Remove "S. STR." string
wood_trait_df$Taxa <- str_replace_all(wood_trait_df$Taxa, "S. STR. ", "")

# Remove "(UNA-DE-GATO, WRIGHT CATCLAW) Synonym: " string
wood_trait_df <- wood_trait_df %>%
  mutate(
    Taxa = str_replace(
      Taxa,
      fixed("(UNA-DE-GATO, WRIGHT CATCLAW)|"),
      ""
    ) |> 
      str_squish()
  )

# Replace "x _" with "x__" to retain crosses
wood_trait_df$Taxa <- str_replace_all(wood_trait_df$Taxa, " x ", " x_")

# Create function to construct columns for family (first word), genus (first word not in all caps) and species from initial taxa column
extract_taxa <- function(taxa) {
  # normalize: replace | with space, compress whitespace
  toks <- str_split(str_squish(str_replace_all(taxa, "\\|", " ")), "\\s+")[[1]]
  if (length(toks) == 0) {
    return(tibble(Family = NA_character_, Genus = NA_character_, Species = NA_character_))
  }
  
  # index of first ALL-CAPS token (the Family)
  fam_idx <- which(str_detect(toks, "^[A-Z]+$"))[1]
  if (is.na(fam_idx)) {
    return(tibble(Family = NA_character_, Genus = NA_character_, Species = NA_character_))
  }
  
  Family <- toks[fam_idx]
  
  # advance through any contiguous ALL-CAPS tokens after family (skip subfamily, tribe, etc.)
  last_allcaps <- fam_idx
  while (last_allcaps + 1 <= length(toks) && str_detect(toks[last_allcaps + 1], "^[A-Z]+$")) {
    last_allcaps <- last_allcaps + 1
  }
  
  # genus is the first token after that block
  genus_idx <- last_allcaps + 1
  if (genus_idx > length(toks)) {
    return(tibble(Family = Family, Genus = NA_character_, Species = NA_character_))
  }
  Genus <- toks[genus_idx]
  
  # species = "Genus <next_token>" if next token exists, else just Genus
  if (genus_idx + 1 <= length(toks)) {
    Species <- paste(Genus, toks[genus_idx + 1])
  } else {
    Species <- Genus
  }
  
  tibble(Family = Family, Genus = Genus, Species = Species)
}

# Apply to wood traits data frame
wood_trait_taxa_df <- wood_trait_df %>%
  mutate(Taxa = str_squish(str_replace_all(Taxa, "\\|", " "))) %>%    # normalize
  bind_cols(map_dfr(.$Taxa, extract_taxa))

# Remove non-specified species - rows containing in species column "SPP.", "SP.", "group", "sect." or "subsect."
wood_trait_taxa_df <- wood_trait_taxa_df %>%
  filter(
    !str_detect(
      Species,
      regex("spp\\.|sp\\.|group|sect\\.|subsect\\.", ignore_case = TRUE)
    )
  )

# Make all empty cells NA
wood_trait_taxa_df <- wood_trait_taxa_df %>%
  mutate(across(where(is.character), ~ na_if(.x, "")))

# Trees present in both Europe and North American datasets will be repeated

# Add indices for filtering criteria
wood_trait_taxa_df <- mutate(wood_trait_taxa_df,
                             Index = 1:nrow(wood_trait_taxa_df))

# Run through duplicate species and select record that is most detailed, with reason
# Prioritise those with most detailed analyses first, then confidence of analysis
duplicate_summary <- wood_trait_taxa_df %>%
  group_by(dataset, Species) %>%
  filter(n() > 1) %>%
  summarize(
    count = n(),
    indices = list(Index),
    sample_rows = list(paste0(Index, ": ", Species)[1:min(5,n())]), # small sample to help review
    .groups = "drop"
  ) %>%
  arrange(Species, dataset)

# View the summary to review duplicates
duplicate_summary$Species

# [1] "Acer maximowiczianum"  remove 831, less detailed     
# [2] "Alangium chinense" remove 222, less detailed         
# [3] "Alnus incana" remove 114, less detailed (missing micrographs)              
# [4] "Balanites aegyptiaca"  remove 987, less detailed     
# [5] "Berchemia floribunda" remove 653, less detailed     
# [6] "Berchemiella berchemiifolia" remove 655, less detailed
# [7] "Calycanthus floridus" remove 1053, 1055, 1053 less detailed, and 1055 very similar      
# [8] "Celtis laevigata" remove 1058, less detailed          
# [9] "Chrysojasminum fruticans" remove 576, less detailed  
# [10] "Citrus x_aurantium" remove 770, less detailed        
# [11] "Cornus sericea" remove 1075, less confident            
# [12] "Cydonia oblonga" remove 688, less detailed           
# [13] "Dendropanax trifidus" remove 84, 85, less detailed      
# [14] "Diospyros ferrea" remove 238, less detailed          
# [15] "Ehretia acuminata"  remove 149, less detailed        
# [16] "Hedera helix"  remove 90, less detailed             
# [17] "Juglans regia" remove 371, less detailed             
# [18] "Kalopanax septemlobus" remove 93, less detailed     
# [19] "Koelreuteria bipinnata" remove 859, less detailed    
# [20] "Laurus nobilis" remove 410, less detailed            
# [21] "Manilkara zapota" remove 1260, less detailed          
# [22] "Melia azedarach" remove 531, less confident           
# [23] "Morus alba" remove 555, less detailed                
# [24] "Nerium oleander" remove 43, less detailed           
# [25] "Olea europaea" remove 599, less confident             
# [26] "Pistacia terebinthus" remove 18, less detailed      
# [27] "Ptelea trifoliata" remove 1229, less detailed         
# [28] "Pterocarya macroptera" remove 375, less confident     
# [29] "Ricinus communis" remove 306, less confident          
# [30] "Sambucus nigra" almost identical, remove 1291            
# [31] "Sambucus racemosa" remove 1292, less detailed         
# [32] "Symplocos anomala"  remove 895, less detailed        
# [33] "Symplocos sumuntia" remove 912, less detailed        
# [34] "Tamarix nilotica"  remove 920, less detailed         
# [35] "Tetradium glabrifolium" remove 781, less detailed    
# [36] "Ulmus davidiana" remove 954, less detailed           
# [37] "Ungnadia speciosa" remove 1258, less detailed         
# [38] "Viburnum odoratissimum" remove 974, less detailed    
# [39] "Viburnum opulus" remove 5, less detailed           
# [40] "Vitex agnus-castus" remove 390, less detailed       
# [41] "Vitis vinifera" remove 984, less detailed            
# [42] "Zelkova serrata" remove 964, less detailed           
# [43] "Ziziphus jujuba" remove 670, less detailed

# Choose duplicates to remove
duplicates_to_remove <- c(831, 222, 987, 653, 655, 1053, 1055, 1058,
                          576, 770, 1075, 688, 84, 85, 238, 149, 371,
                          93, 859, 410, 1260, 531, 555, 599, 18, 1229,
                          375, 306, 1291, 895, 912, 920, 781, 954,
                          1258, 974, 5, 390, 984, 964, 670)

# Remove the selected duplicate records by their Index values.
wood_traits_nodup_df <- wood_trait_taxa_df %>%
  filter(!Index %in% duplicates_to_remove)

# The GeoTIFF stores the climate class number as its cell value and its RGB
# colour in a colour table. Match that RGB colour to the imported legend.
raster_coltab_df <- as.data.frame(terra::coltab(kg_raster)[[1]]) %>%
  transmute(
    value = as.integer(value),
    rgb = paste(red, green, blue)
  )

koppen_geiger_lookup <- raster_coltab_df %>%
  inner_join(koppen_geiger_legend_df, by = c("value", "rgb")) %>%
  select(value, koppen_geiger)

host_pref_df <- host_pref_df %>%
  mutate(
    .climate_value = as.integer(terra::values(kg_raster)[
      terra::cellFromXY(kg_raster, cbind(lon, lat)),
      1
    ])
  ) %>%
  left_join(koppen_geiger_lookup, by = c(".climate_value" = "value")) %>%
  select(-.climate_value)

# Check duplicates removed
duplicate_check <- wood_traits_nodup_df %>%
  group_by(Species) %>%
  filter(n() > 1) %>%
  summarize(
    count = n(),
    indices = list(Index),
    sample_rows = list(paste0(Index, ": ", Species)[1:min(5,n())]), # small sample to help review
    .groups = "drop"
  ) 

# Query GBIF for species distributions ------------------------------------
# Köppen-Geiger classes represented in the host-preference studies.
study_kg_classes <- host_pref_df %>%
  pull(koppen_geiger) %>%
  str_split(";", simplify = FALSE) %>%
  unlist() %>%
  str_squish() %>%
  discard(~ is.na(.x) || .x == "") %>%
  unique()

# Approximate geographic limits for the GBIF searches, expressed as
# c(west, south, east, north). The European limits span Ireland to Iran and
# the Mediterranean to Scandinavia, as specified for this analysis.
europe_gbif_bbox <- c(-10.7, 25.0, 63.3, 71.5)
california_gbif_bbox <- c(-124.5, 32.5, -114.1, 42.1)

# Cache one result per dataset-species query so completed queries are not
# repeated on subsequent runs.
gbif_cache_dir <- "gbif_cache"
if (!dir.exists(gbif_cache_dir)) dir.create(gbif_cache_dir)

gbif_cache_file <- function(species, dataset) {
  cache_name <- paste(dataset, species, sep = "__") %>%
    str_replace_all("[^A-Za-z0-9]+", "_")
  file.path(gbif_cache_dir, paste0(cache_name, ".rds"))
}

get_kg_classes_from_occurrences <- function(occurrences,
                                            target_classes = study_kg_classes) {
  if (nrow(occurrences) == 0 || length(target_classes) == 0) {
    return(character())
  }

  occurrence_points <- occurrences %>%
    transmute(
      lon = as.numeric(decimalLongitude),
      lat = as.numeric(decimalLatitude)
    ) %>%
    filter(is.finite(lon), is.finite(lat)) %>%
    distinct()

  if (nrow(occurrence_points) == 0) return(character())

  raster_values <- terra::extract(
    kg_raster,
    occurrence_points[, c("lon", "lat")]
  )[[2]]

  koppen_geiger_lookup %>%
    filter(
      value %in% raster_values,
      koppen_geiger %in% target_classes
    ) %>%
    distinct(koppen_geiger) %>%
    pull(koppen_geiger)
}

# Query GBIF in pages so that all available coordinate records within the
# selected geographic limits can be considered, up to a maximum of 5,000
# records. Pagination stops early once all study classes have been found.
get_gbif_occurrences <- function(species, dataset, page_size = 300,
                                 max_records = 5000,
                                 target_classes = study_kg_classes) {
  gbif_species <- str_replace(species, " x_", " x ")

  query_args <- list(
    scientificName = gbif_species,
    hasCoordinate = TRUE,
    fields = c("decimalLongitude", "decimalLatitude"),
    limit = page_size
  )

  if (dataset == "European") {
    query_args$geometry <- europe_gbif_bbox
  } else if (dataset == "North American") {
    query_args$stateProvince <- "California"
    query_args$geometry <- california_gbif_bbox
  } else {
    return(tibble())
  }

  all_records <- list()
  observed_classes <- character()
  n_records <- 0
  start <- 0

  repeat {
    records_remaining <- max_records - n_records
    if (records_remaining <= 0) break
    query_args$limit <- min(page_size, records_remaining)

    result <- tryCatch(
      do.call(rgbif::occ_search, c(query_args, start = start)),
      error = function(error) {
        warning(
          "GBIF query failed for ", species, " (", dataset, "): ",
          conditionMessage(error)
        )
        NULL
      }
    )

    if (is.null(result) || is.null(result$data) || nrow(result$data) == 0) {
      break
    }

    all_records[[length(all_records) + 1]] <- result$data
    n_records <- sum(vapply(all_records, nrow, integer(1)))
    observed_classes <- union(
      observed_classes,
      get_kg_classes_from_occurrences(result$data, target_classes)
    )

    if (isTRUE(result$meta$endOfRecords) ||
        nrow(result$data) < query_args$limit ||
        n_records >= max_records ||
        (length(target_classes) > 0 &&
         all(target_classes %in% observed_classes))) {
      break
    }

    start <- start + query_args$limit
  }

  if (length(all_records) == 0) tibble() else bind_rows(all_records)
}

# Convert GBIF coordinates to Köppen-Geiger classes and retain only classes
# represented by the host-preference studies.
get_species_kg_classes <- function(species, dataset) {
  occurrences <- get_gbif_occurrences(species, dataset)

  if (nrow(occurrences) == 0) {
    return(tibble(
      dataset = dataset,
      Species = species,
      koppen_geiger = NA_character_
    ))
  }

  observed_classes <- get_kg_classes_from_occurrences(occurrences) %>%
    sort()

  tibble(
    dataset = dataset,
    Species = species,
    koppen_geiger = if (length(observed_classes) == 0) {
      NA_character_
    } else {
      paste(observed_classes, collapse = "; ")
    }
  )
}

species_to_query <- wood_traits_nodup_df %>%
  distinct(dataset, Species) %>%
  filter(
    dataset %in% c("European", "North American"),
    !is.na(Species),
    Species != ""
  )

gbif_progress <- utils::txtProgressBar(
  min = 0,
  max = nrow(species_to_query),
  style = 3
)

gbif_kg_results <- vector("list", nrow(species_to_query))

for (i in seq_len(nrow(species_to_query))) {
  species <- species_to_query$Species[i]
  dataset <- species_to_query$dataset[i]
  cache_file <- gbif_cache_file(species, dataset)

  if (file.exists(cache_file)) {
    gbif_kg_results[[i]] <- readRDS(cache_file)
  } else {
    gbif_kg_results[[i]] <- get_species_kg_classes(species, dataset)
    saveRDS(gbif_kg_results[[i]], cache_file)
  }

  utils::setTxtProgressBar(gbif_progress, i)
}

close(gbif_progress)

gbif_kg_results <- bind_rows(gbif_kg_results)

write.csv(gbif_kg_results, "gbif_kg_results.csv")

# GO FROM HERE
wood_traits_nodup_df <- wood_traits_nodup_df %>%
  select(-any_of("koppen_geiger")) %>%
  left_join(gbif_kg_results, by = c("dataset", "Species"))


# Calculate wood traits ---------------------------------------------------
# Rename porosity columns
wood_traits_nodup_df <- wood_traits_nodup_df %>%
  rename(ring_porous = `3 - Wood ring-porous`,
         semi_porous = `4 - Wood semi-ring-porous`,
         diffuse_porous = `5 - Wood diffuse-porous`)

# Calculate porosity score for each species, weighted by variability (0.5)
porosity_df <- wood_traits_nodup_df %>%
  pivot_longer(cols = c(ring_porous, semi_porous, diffuse_porous),
               names_to = "type",
               values_to = "val") %>%
  mutate(
    # porosity score per type
    score = case_when(
      type == "ring_porous"   ~ 0,
      type == "semi_porous"   ~ 0.5,
      type == "diffuse_porous"~ 1
    ),
    # weight: 0 if missing, 0.5 if contains "v", otherwise 1
    weight = case_when(
      is.na(val)            ~ 0,
      grepl("v", val, fixed = TRUE) ~ 0.5,
      TRUE                  ~ 1
    )
  ) %>%
  group_by(Species) %>%
  summarise(
    wood_porosity = if (sum(weight) == 0) NA_real_ else sum(score * weight) / sum(weight),
    .groups = "drop"
  )

# Rejoin to original dataframe
wood_traits_porosity_df <- full_join(porosity_df, wood_traits_nodup_df)

# Plot distribution of porosity across species
ggplot(data = wood_traits_porosity_df, aes(x = wood_porosity)) +
  geom_histogram() +
  labs(x = "Wood porosity score", y = "Number of species") +
  theme_bw()

# Calculate for each genus - run through where not in agreement across genus? or calculate mean across genus?
porosity_genus_df <- wood_traits_porosity_df %>% 
  group_by(Genus) %>%
  summarise(
    mean_porosity = mean(wood_porosity, na.rm = TRUE),
    sd_porosity   = sd(wood_porosity, na.rm = TRUE),
    n = n(), # Sample size
    .groups = "drop"
  )

# Plot distribution of porosity across genera
ggplot(data = porosity_genus_df, aes(x = mean_porosity)) +
  geom_histogram() +
  labs(x = "Mean wood porosity score", y = "Number of genera") +
  theme_bw()
