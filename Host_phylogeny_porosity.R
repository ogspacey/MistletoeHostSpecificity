### Host phylogeny and porosity
### Oliver G. Spacey
### 2025.12.04

# Code to plot a phylogeny of hosts of Viscum album subsp. album across angiosperms
# And their traits, e.g., wood porosity

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load required packages
library(tidyverse) # Data wrangling
library(rotl)      # Phylogenetic analysis
library(ape)       # Phylogenetic analysis
library(phytools)  # Phylogenetic analysis
library(phylolm)   # Phylogenetic analysis
library(ggplot2)   # Data visualisation


# Load data ---------------------------------------------------------------
# Load wood trait datasets
# Load european trees
european_df <- read.csv("European_IW_data.csv", header = TRUE)
# Load North American trees
n_america_df <- read.csv("N_America_IW_data.csv", header = TRUE)
# Make all columns characters before combining
for(i in 1:ncol(european_df)){
  european_df[, i] <- as.character(european_df[, i])
}
for(i in 1:ncol(n_america_df)){
  n_america_df[, i] <- as.character(n_america_df[, i])
}
  
# Combine wood trait datasets
wood_traits_df <- bind_rows(european_df, n_america_df)

# Load hosts of Viscum album subsp album
host_genera_df <- read.csv("viscum_album_host_genera.csv", header = TRUE) %>%
                  mutate(host = 1) # Add column to designate these are host genera

# Load preference data - LATER

# Wrangle data ------------------------------------------------------------
# Find obscure strings and remove
# Remove "S. STR." string
wood_traits_df$Taxa <- str_replace_all(wood_traits_df$Taxa, "S. STR. ", "")

# Remove "(UNA-DE-GATO, WRIGHT CATCLAW) Synonym: " string
wood_traits_df <- wood_traits_df %>%
  mutate(
    Taxa = str_replace(
      Taxa,
      fixed("(UNA-DE-GATO, WRIGHT CATCLAW)|"),
      ""
    ) |> 
      str_squish()
  )

# Replace "x _" with "x__" to retain crosses
wood_traits_df$Taxa <- str_replace_all(wood_traits_df$Taxa, " x ", " x_")

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
wood_traits_taxa_df <- wood_traits_df %>%
  mutate(Taxa = str_squish(str_replace_all(Taxa, "\\|", " "))) %>%    # normalize
  bind_cols(map_dfr(.$Taxa, extract_taxa))

# Remove non-specified species - rows containing in species column "SPP.", "SP.", "group", "sect." or "subsect."
wood_traits_taxa_df <- wood_traits_taxa_df %>%
  filter(
    !str_detect(
      Species,
      regex("spp\\.|sp\\.|group|sect\\.|subsect\\.", ignore_case = TRUE)
    )
  )

# Trees present in both Europe and North American datasets will be repeated - remove
# Make all empty cells NA
wood_traits_taxa_df <- wood_traits_taxa_df %>%
  mutate(across(everything(), ~ na_if(.x, "")))

# Remove columns specifying geographical location (can retrieve at later time from original datasets)
wood_traits_taxa_df <- select(wood_traits_taxa_df, -c(159:183))

# Keep only distinct rows
wood_traits_taxa_df <- wood_traits_taxa_df %>% distinct()

# Add indices for filtering criteria
wood_traits_taxa_df <- mutate(wood_traits_taxa_df,
                              Index = 1:nrow(wood_traits_taxa_df))

# Run through duplicate species and select record that is most detailed, with reason
# Prioritise those with most detailed analyses first, then confidence of analysis
duplicate_summary <- wood_traits_taxa_df %>%
  group_by(Species) %>%
  filter(n() > 1) %>%
  summarize(
    count = n(),
    indices = list(Index),
    sample_rows = list(paste0(Index, ": ", Species)[1:min(5,n())]), # small sample to help review
    .groups = "drop"
  ) 

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
duplicates_removed <- c(831, 222, 114, 987, 653, 655, 1053, 1055, 1058,
                        576, 770, 1075, 688, 84, 85, 238, 149, 90, 371,
                        93, 859, 410, 1260, 531, 555, 43, 599, 18, 1229,
                        375, 306, 1291, 1292, 895, 912, 920, 781, 954, 
                        1258, 974, 5, 390, 984, 964, 670)
wood_traits_nodup_df <- wood_traits_taxa_df %>%
                        filter(!Index %in% duplicates_removed)

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

# Rename porosity columns
wood_traits_nodup_df <- wood_traits_nodup_df %>%
                        rename(ring_porous = X3...Wood.ring.porous,
                               semi_porous = X4...Wood.semi.ring.porous,
                               diffuse_porous = X5...Wood.diffuse.porous)

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

# Combine with host genera dataframe
host_porosity_df <- full_join(porosity_genus_df, host_genera_df)

# Remove genera which have no porosity data
host_porosity_df <- host_porosity_df %>%
                    filter(!is.na(mean_porosity))

# Designate all non-host genera as non-hosts
host_porosity_df <- host_porosity_df %>%
  mutate(host = if_else(is.na(host), 0, host))

# Plot distribution of porosity and hosts
ggplot(data = host_porosity_df, aes(x = mean_porosity, fill = as.factor(host))) +
  geom_histogram(position = "identity") +
  labs(x = "Wood porosity score", y = "Number of genera") +
  theme_bw()

# Calculate mean and standard deviation of wood porosity score for hosts and non-hosts
host_summary <- host_porosity_df %>%
  group_by(host) %>%
  summarise(
    mean_porosity = mean(mean_porosity, na.rm = TRUE),
    sd_porosity   = sd(mean_porosity, na.rm = TRUE),
    n = n(),  
    .groups = "drop"
  )

# Perform Wilcoxon test to compare hosts and non-hosts
wilcox.test(
  mean_porosity ~ host,
  data = host_porosity_df,
  exact = FALSE     # recommended if sample sizes are > ~50 or ties present
)

# Plot distribution with boxplot
ggplot(host_porosity_df, aes(x = as.factor(host),
                             y = mean_porosity,
                             fill = as.factor(host))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +    # boxplot
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) + # jittered points
  scale_x_discrete(labels = c("0" = "Non-host", "1" = "Host")) +
  labs(
    x = "",
    y = "Mean porosity (genus level)",
    title = "Genus-level wood porosity for host vs non-host genera"
  ) +
  theme_bw()

# Examine native and inoculated hosts

# Phylogenetic analysis ---------------------------------------------------
# Get phylogeny from Open Tree of Life with genera from subsetted dataset
# Get OTT IDs
taxa <- tnrs_match_names(unique(tolower(host_porosity_df$Genus)), context = "Land plants")

# Map names in OTT to our names
taxon_map <- structure(taxa$search_string, names = taxa$unique_name)

# Get the tree for these taxa
tree <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))])

# Plot tree
plot(tree, show.tip.label = FALSE)

# Replace tip labels with our labels
otl_tips <- strip_ott_ids(tree$tip.label, remove_underscores = TRUE)
tree$tip.label <- taxon_map[ otl_tips ]
tree$node.label <- NULL

# Plot tree with new labels
plot(tree, show.tip.label = TRUE)

# Make genera lower case
host_porosity_df$Genus <- tolower(host_porosity_df$Genus)

# Keep tip labels present in data
keep <- intersect(tree$tip.label, host_porosity_df$Genus)

# Prune tree 
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, keep))
plot(tree_pruned, show.tip.label = FALSE)

# Prune data so matching
pruned_df <- host_porosity_df %>% filter(Genus %in% tree_pruned$tip.label)

# Reorder rows to match tree tip order and name rows accordingly
pruned_df <- pruned_df[match(tree_pruned$tip.label, pruned_df$Genus), ]
rownames(pruned_df) <- pruned_df$Genus

# Add branch lengths using Grafen's method
tree_bl <- compute.brlen(tree_pruned, method = "Grafen")

# Calculate phylogenetic inertia (Pagel's λ) for porosity
pagels_lambda <- phylosig(tree_bl, x = pruned_df$mean_porosity, method = "lambda", test = TRUE)
pagels_lambda

# Model relationship between porosity score and host (phylogenetic logistic regression)
phylo_logit <- phyloglm(host ~ mean_porosity,
                        phy = tree_bl,
                        data = pruned_df,
                        method = "logistic_MPLE")

summary(phylo_logit)

# Plot phylogeny of host genera, their porosity score and host status
# Create vector for porosity values
porosity_vec <- pruned_df$mean_porosity[match(tree_bl$tip.label, pruned_df$Genus)]
names(porosity_vec) <- tree_bl$tip.label

# Create vector for host values
stopifnot(all(tree_bl$tip.label == pruned_df$Genus))
host_vec <- pruned_df$host[match(tree_bl$tip.label, pruned_df$Genus)]

# Plot porosity as continuous trait on phylogeny
cont <- contMap(type = "phylogram", tree_bl, porosity_vec, fsize = c(1, 1))
plot(cont, outline = FALSE, legend = 0.5, cex = 1)

# Plot host tips
host_tips <- which(host_vec == 1)
tiplabels(pch = 21, cex = 1.2, tip = host_tips, offset = 0.03, bg = "red")
