# ==============================================================================================================
#                                               --- START-UP ---
# ==============================================================================================================


# Install packages for spatial data mapping and handling
install.packages(c(
  "dplyr",
  "ggplot2",
  "sysfonts",
  "showtext",
  "here",
  "terra", 
  "geodata", 
  "rnaturalearth", 
  "rnaturalearthdata",
  "dismo", 
  "tidyverse", 
  "rgbif", 
  "raster",
  "caret"
))

# Load packages
library(dplyr)
library(ggplot2)
library(here)
library(terra)
library(geodata)
library(rnaturalearth)
library(dismo)
library(tidyr)
library(rgbif)
library(sysfonts)
library(showtext)
library(caret)

# Create folders 
folders <- c("data/raw", "data/processed", "outputs/maps")
for (f in folders) {
  dir.create(here(f), recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================================================
#                                           --- REGION SET-UP ---
# ==============================================================================================================


# Region presets with lon and lat bounds - add more if desired
REGION_PRESETS <- list(
  "Africa"        = c(-20, 55, -40, 37.5),
  "Europe"        = c(-25, 60, 34, 72),
  "South America" = c(-95, -30, -60, 15),
  "North America" = c(-170, -50, 5, 85),
  "Oceania"       = c(110, 180, -50, 0),
  "Asia"          = c(25, 180, -10, 85),
  "Global"        = c(-180, 180, -90, 90)
)


# ==============================================================================================================
#                                              --- TASK 1 ---
# ==============================================================================================================

# This section involves running GLMs to  predict the present-day distribution of Species 1
# and Species 2 using climate variables. Then, these are used to present maps of their
# current distribution.
# Here, a MAIN function plots an SDM for the user's Species 1 and Species 2 in a user-defined region
# (Africa, Europe, South America, North America, Oceania, Asia, Global). The user can pick what 
# bioclimatic variables are used, or let it be decided algorithmically.

# 4 HELPER functions have been made to split the generation of an SDM into 4 parts:
# -- 1: Download and clean occurrence, ocean, and climate data.
# -- 2 (OPTIONAL): Pick bioclimatic variables based on multicollinearity and stepwise selection.
# -- 3: Prepare presence and background/pseudoabsence points, train/test data, fit GLM, and evaluate.
# -- 4: Predict the habitat suitability for both species separately across the region and map this.


# HELPER FUNCTION 1 - DOWNLOAD & CLEAN DATA
# ==============================================================================================================
get_species_data <- function(sp_name, region, region_bounds, ocean_vect, bioclim_global) {
  
  sp_filename <- gsub(" ", "_", sp_name)
  
  # -1- Define file path and download occurrence data ----------------------------------------------------------
  sp_file <- here("data", "raw", paste0(sp_filename, ".rds"))
  
  if (!file.exists(sp_file)) {
    message("Downloading occurrences from GBIF (may take a few minutes)...")
    # Limit download to 10,000 data points
    occ <- occ_search(scientificName = sp_name, hasCoordinate = TRUE, limit = 10000)
    occ_df <- occ$data
    saveRDS(occ_df, sp_file)
  } else {
    # Load if already exists
    occ_df <- readRDS(sp_file)
  }
  
  # [SAFETY CHECK] Checks if there are any records at all
  if (is.null(occ_df) || nrow(occ_df) == 0) {
    message(paste("SKIP: No records found for", sp_name, "- check spelling or internet connection!"))
    return(NULL)
  }
  
  
  # -1.1- Extract and clean coordinates in chosen region ---
  coords_region <- occ_df %>% 
    dplyr::select(decimalLongitude, decimalLatitude) %>%
    rename(lon = decimalLongitude, lat = decimalLatitude) %>%
    na.omit() %>%
    filter(lon >= region_bounds[1], lon <= region_bounds[2],
           lat >= region_bounds[3], lat <= region_bounds[4])
  
  # [SAFETY CHECK] Check raw occurrence in region before any heavy processing
  if (nrow(coords_region) < 50) {
    message(paste("SKIP: Not enough raw points for", sp_name, "in", region,
                  "(Found:", nrow(coords_region), "- Minimum required: 50)"))
    return(NULL)
  }
  cat("Records in ", region, ":", nrow(coords_region), "\n")
  
  
  # -1.3- Remove occurrence points in the ocean ---
  
  # Convert coordinates to SpatVector
  species_vect <- vect(coords_region, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Make sure CRS matches ocean
  crs(species_vect) <- crs(ocean_vect)
  
  # Check every point against every ocean polygon
  ocean_intersects <- relate(species_vect, ocean_vect, relation = "intersects")
  
  # Create a logical vector that returns TRUE if a point is in any ocean polygon
  is_ocean <- apply(ocean_intersects, 1, any)
  
  # Store the points that DO NOT intersect with the ocean polygons
  species_land_vect <- species_vect[!is_ocean, ]
  
  
  # -2- Environmental data and WorldClim bioclimatic variables -------------------------------------------------
  
  # -2.1- Define a study extent and reduce sampling bias ---
  
  # Create extent from the bounds passed into the function
  region_ext <- ext(region_bounds[1], region_bounds[2], region_bounds[3], region_bounds[4])
  
  # Crop the global bioclim (passed as argument) to this extent
  bioclim_crop <- crop(bioclim_global, region_ext)
  
  # Identify the cell number for each point
  cells <- cellFromXY(bioclim_crop, geom(species_land_vect)[, c("x", "y")])
  
  # Remove duplicate points
  species_land_vect <- species_land_vect[!duplicated(cells), ]
  
  cat("Points remaining after spatial thinning:", length(species_land_vect), "\n")
  
  
  # -2.2- Extract raster values at occurrence points ---
  
  # Ensure point CRS matches rasters
  species_pts <- project(species_land_vect, crs(bioclim_crop))
  
  # Extract climate values at points
  clim_vals <- terra::extract(bioclim_crop, species_pts)
  
  # Combine coordinates with climate values
  species_data <- bind_cols(
    as_tibble(crds(species_pts)) %>% rename(lon = x, lat = y),
    as_tibble(clim_vals)[,-1]
  ) %>% drop_na()
  
  # [SAFETY CHECK] Check processed occurrence points in region
  if (nrow(species_data) < 50) {
    message(paste("SKIP: Not enough processed occurrence points for", sp_name,
                  "in", region, ". Only", nrow(species_data), "valid points found."))
    return(NULL)
  }
  cat("Final number of dataset rows:", nrow(species_data), "\n")
  
  # Save processed data
  write.csv(species_data, here("data", "processed", paste0(sp_filename, ".csv")))
  
  # Return data and cropped climate layers
  return(list(
    data = species_data, 
    clim_crop = bioclim_crop))
}
# ==============================================================================================================




# HELPER FUNCTION 2 - BIOCLIMATIC VARIABLE SELECTION
# ==============================================================================================================
# Function which takes a presence/background points + bioclim vars 1-19 to find most parsimonious mode
# The user could run this function, and then pick what variables make most sense based on previous
# knowledge.
bioclim_selection <- function(input_data) {
  
  # Only look at 'bio' columns for correlation
  bio_cols <- grep("bio", names(input_data), value = TRUE)
  
  
  # -1- Checking for multicollinearity -------------------------------------------------------------------------
  
  # Calculate correlation matrix
  cor_matrix <- cor(input_data[, bio_cols], method = "spearman")
  
  # Find the most highly correlated variables - these will be removed
  high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.7)
  
  # Store names of weakly correlated variables to keep
  if (length(high_cor_vars) > 0) {
    clean_vars <- bio_cols[-high_cor_vars]
  } else {
    clean_vars <- bio_cols
  }
  
  # -2- Stepwise BIC selection ---------------------------------------------------------------------------------
  
  # Fit a model using only weakly correlated variables
  form_start <- as.formula(paste("presence ~", paste(clean_vars, collapse = "+")))
  
  # Suppress warnings for initial fit
  full_model <- suppressWarnings(glm(form_start, data = input_data, family = binomial))
  
  # Run stepwise selection which uses BIC to compare different subsets of the model
  # BIC is stricter than AIC, particularly for bigger sample size (n)
  n <- nrow(input_data)
  best_model <- step(full_model, direction = "both", trace = 0, k = log(n))
  
  # Extract the names of the 'winning' variables (and remove the intercept)
  best_vars <- names(coef(best_model))[-1]
  
  return(best_vars)
  
}
# ==============================================================================================================




# HELPER FUNCTION 3 - PRESENCE/BACKGROUND POINTS, FIT GLM & EVALUATE
# ==============================================================================================================
fit_eval_glm <- function(species_data, bioclim_crop, user_predictors = NULL) {
  
  # -1- Generate presence/pseudoabsence points -----------------------------------------------------------------
  
  # Dynamic sample size for background points (either 1000 or number of presences if > 1000)
  bg_n <- max(1000, nrow(species_data))
  
  # Sample random points from the first layer of the cropped climate data
  set.seed(123)
  bg_pts <- spatSample(bioclim_crop[[1]], size = bg_n, method = "random", 
                       na.rm = TRUE, as.points = TRUE, values = FALSE)
  
  # Extract climate values for these background points
  bg_clim <- terra::extract(bioclim_crop, bg_pts)
  
  # Format background point data
  bg_coords <- as.data.frame(crds(bg_pts))
  colnames(bg_coords) <- c("lon", "lat")
  
  background_data <- bind_cols(bg_coords, as_tibble(bg_clim)[,-1]) %>%
    mutate(presence = 0) %>%
    drop_na()
  
  cat("Background points after NAs removed:", nrow(background_data), "\n")
  
  # Format presence data
  presence_data <- species_data %>% mutate(presence = 1)
  
  
  # -1.1- Split data into training (70%) and testing (30%) data ---
  
  # Split presence data
  k = 0.7
  pres_idx <- sample(nrow(presence_data), size = floor(k * nrow(presence_data)))
  train_pres <- presence_data[pres_idx, ]
  test_pres <- presence_data[-pres_idx, ]
  
  # Split background data
  bg_idx <- sample(nrow(background_data), size = floor(k * nrow(background_data)))
  train_bg <- background_data[bg_idx, ]
  test_bg <- background_data[-bg_idx, ]
  
  # Combine for training
  train_data <- bind_rows(train_pres, train_bg)
  
  
  # -1.2- Fit GLM ---
  
  current_vars <- user_predictors
  
  # If no variables provided, run bioclim_selection() function
  if (is.null(current_vars)) {
    message("Calculating most parsimonious model...")
    # Call function to check for multicollinearity and perform stepwise AIC selection
    current_vars <- bioclim_selection(train_data)
  }
  
  # [SAFETY CHECK] Check if any bioclim vars were actually stored
  if (length(current_vars) == 0) {
    stop("Variable selection failed. No predictors found.")
  }
  
  message(paste("Fitting GLM with variables: ", paste(current_vars, collapse = ", ")))
  
  # Construct formula dynamically
  formula_str <- paste("presence ~ ", paste(current_vars, collapse = " + "))
  model_formula <- as.formula(formula_str)
  
  # Fit model using training data
  sdm_model <- glm(model_formula, data = train_data, family = binomial)
  
  
  # -1.4- Evaluate model ---
  
  # Evaluate using held-out test data
  eval_res <- evaluate(
    p = test_pres %>% dplyr::select(all_of(current_vars)),
    a = test_bg %>% dplyr::select(all_of(current_vars)),
    model = sdm_model
  )
  cat("AUC score:", eval_res@auc, "\n")
  
  # Return model, bioclim variables, AUC score
  return(list(
    model = sdm_model, 
    vars = current_vars, 
    auc = eval_res@auc,
    train_data = train_data,
    eval_res = eval_res))
}
# ==============================================================================================================




# HELPER FUNCTION 4 - PREDICTION & MAPPING SDM
# ==============================================================================================================
predict_and_map <- function(sdm_model, bioclim_crop, current_vars, sp_name, species_data, eval_res) {
  
  # -1- Predict and Map ----------------------------------------------------------------------------------------
  
  # Ensure raster has those layers
  if (!all(current_vars %in% names(bioclim_crop))) {
    stop("ERROR: Predictor layers not found in bioclim_crop.")
  }
  
  # Predict probability of occurrence across the full study extent
  prediction <- terra::predict(bioclim_crop[[current_vars]], sdm_model, type = "response")
  names(prediction) <- "suitability"
  
  # Plot continuous suitability map
  plot(prediction, main = paste("Predicted suitability:", sp_name))
  
  # Re-create vector from data just for plotting points
  # (We do this because we can't easily pass the vector object from Helper 1 to Helper 4)
  species_vect_plot <- vect(species_data, geom=c("lon", "lat"), crs=crs(prediction))
  points(species_vect_plot, pch = 16, cex = 0.5, col = "black")
  
  # Calculate threshold using evaluation results so that probability -> presence/absence
  thr <- dismo::threshold(eval_res, 'prevalence')
  
  # Create binary Presence/Absence map
  # This will be needed for later tasks to calculate overlap.
  prediction_pa <- prediction > thr
  names(prediction_pa) <- "presence_absence"
  
  message("Completed generation of current SDMs for ", sp_name)
  cat("-----------------------------------------------------------------------\n")
  
  return(list(
    continuous_map = prediction,
    binary_map = prediction_pa,
    threshold = thr))
}
# ==============================================================================================================




# MAIN TASK 1 FUNCTION - CURRENT SDM GENERATION
# ==============================================================================================================
# Function generates current SDMs for a pair of species in a user-defined region
run_current_sdm <- function(species1, species2, region, predictor_list) {
  
  # -1- Load static data (Ocean & Climate) ---------------------------------------------------------------------
  message("Loading environmental and ocean data (may take a few seconds)...")
  
  # -1.1- Prepare Ocean Data ---
  ocean_data_dir <- here("data", "raw", "ocean")
  if (!dir.exists(ocean_data_dir)) dir.create(ocean_data_dir)
  URL <- "https://naturalearth.s3.amazonaws.com/110m_physical/ne_110m_ocean.zip"
  zip_file <- file.path(ocean_data_dir, basename(URL))
  
  if (!file.exists(zip_file)) download.file(URL, zip_file)
  if (length(list.files(ocean_data_dir, pattern = "\\.shp$")) == 0) unzip(zip_file, exdir = ocean_data_dir)
  
  # Read Ocean Shapefile
  shp_file <- list.files(ocean_data_dir, pattern = "\\.shp$", full.names = TRUE)[1]
  ocean <- vect(shp_file)
  
  
  # -1.2- Prepare Climate Data ---
  bioclim_dir <- here("data", "raw", "worldclim")
  if(!dir.exists(bioclim_dir)) dir.create(bioclim_dir)
  
  # Download 19 bioclim variables
  bioclim_global <- worldclim_global(var = "bio", res = 10, path = bioclim_dir)
  names(bioclim_global) <- paste0("bio", 1:19)
  
  
  # -1.3- Get Region Bounds ---
  if (!region %in% names(REGION_PRESETS)) {
    stop("Region not found. Please choose from the preset list.")
  }
  bounds <- REGION_PRESETS[[region]]
  
  
  # -2- Iterate through species using Helper Functions ---------------------------------------------------------
  
  # Prepare lists and dataframes for results output
  current_sdm_maps <- list()
  current_models <- list()
  current_data <- list()
  
  results_stats <- data.frame(
    Species = character(),
    Region = character(),
    N_Points = numeric(),
    Bioclim_Vars = character(),
    AUC = numeric(),
    stringsAsFactors = FALSE
  )
  
  chosen_species <- c(species1, species2)
  
  for (sp_name in chosen_species) {
    
    message(paste("Processing:", sp_name, "in", region))
    
    # -- STEP 1: Download, clean, filter ocean, crop climate, and extract points ---
    prepared_data <- get_species_data(sp_name, region, bounds, ocean, bioclim_global)
    
    # Skip species if there is an error or low sample size
    if (is.null(prepared_data)) next 
    
    # Extract the specific outputs we need for the next steps
    sp_data      <- prepared_data$data
    clim_cropped <- prepared_data$clim_crop
    
    
    # -- STEP 2: Generate background points, split data, select bioclim, fit GLM ---
    glm_results <- fit_eval_glm(sp_data, clim_cropped, predictor_list[[sp_name]])
    
    # Extract outputs
    final_model <- glm_results$model
    final_vars  <- glm_results$vars
    final_auc   <- glm_results$auc
    final_eval  <- glm_results$eval_res
    
    
    # -- STEP 3: Predict the model onto the raster and create the plot ---
    suitability_map <- predict_and_map(final_model, clim_cropped, final_vars, sp_name, sp_data, final_eval)
    
    
    # -- STEP 4: Store Results ---
    current_sdm_maps[[sp_name]] <- suitability_map
    current_models[[sp_name]] <- glm_results$model
    current_data[[sp_name]] <- glm_results$train_data
    
    results_stats <- rbind(results_stats, data.frame(
      Species = sp_name,
      Region = region,
      N_Points = nrow(sp_data),
      Bioclim_Vars = paste(final_vars, collapse = ", "),
      AUC = final_auc
    ))
    
  }
  
  message("All SDMs generated successfully. Returning map list...")
  return(list(
    maps = current_sdm_maps,
    models = current_models,
    data = current_data,
    stats = results_stats))
  
}
# ==============================================================================================================




# ==============================================================================================================
#                                               --- TASK 1 EXECUTION ---
# ==============================================================================================================

# Define species 1 and species 2 - USER CHANGE
sp1 <- "Loranthus europaeus"
sp2 <- "Quercus petraea"

# Define chosen bioclimatic variables for each species - USER CHANGE (IF DESIRED)
sp_predictors <- list()
sp_predictors[[sp1]] <- NULL
sp_predictors[[sp2]] <- NULL

# Results
sdm_results <- run_current_sdm(sp1, sp2, "Europe", sp_predictors)
my_maps <- sdm_results$maps
my_results <- lapply(sdm_results$models, summary)

print(sdm_results$stats)
print(my_results)








# ==============================================================================================================
#                                              --- TASK 2 ---
# ==============================================================================================================

# This section involves plotting the overlap in distribution of Species 1 and Species 2, and
# devising and calculating a metric for the degree of overlap between their ranges.
# Task 2 will be converted into one general function that can take any Species 1 and 2, and the
# results outputed from the previous task.


# MAIN TASK 2 FUNCTION - DISTRIBUTION OVERLAP METRIC AND PLOT
# ==============================================================================================================
calculate_and_map_overlap <- function(species1, species2, sdm_results) {
  
  # -1- Calculating overlap of Species 1 and Species 2 using binary maps ---------------------------------------
  
  # -1.1- Extract binary rasters from Task 1 predict_and_map() function ---
  sp1_bin <- sdm_results$maps[[species1]]$binary_map
  sp2_bin <- sdm_results$maps[[species2]]$binary_map
  
  
  # -1.2- Calculate the degree of overlap (intersection/union)
  
  # Intersection (both species present)
  overlap_rast <- sp1_bin & sp2_bin
  
  # Union (at least one species is present)
  union_rast <- sp1_bin | sp2_bin
  
  # Count pixels using terra::global() - na.rm = TRUE ignores NAs/ocean points
  overlap_pixels <- terra::global(overlap_rast, "sum", na.rm = TRUE)[1, 1]
  union_pixels <- terra::global(overlap_pixels, "sum", na.rm = TRUE)[1, 1]
  
  # Calculate percentage of overlap 
  overlap_metric <- (overlap_pixels / union_pixels) * 100
  
  cat("Degree of overlap of", species1, "and", species2, ":", round(overlap_metrix, 2), "%\n")
  
  
  
}





























































