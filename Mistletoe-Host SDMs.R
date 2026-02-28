# Stephanie Pearce
# Data accessed: 28/02/2026 

# ==============================================================================================================
#                                               --- START-UP ---
# ==============================================================================================================


# -1- Define all required packages for SDM pipeline ---
required_packages <-c(
  "dplyr",
  "ggplot2",
  "sysfonts",
  "showtext",
  "here",
  "terra",
  "tidyterra",
  "geodata", 
  "rnaturalearth", 
  "rnaturalearthdata",
  "dismo", 
  "tidyverse", 
  "rgbif", 
  "raster",
  "caret",
  "sf",
  "ggspatial"
)

# Install all the packages the user is missing
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, dependencies = TRUE)
}

# Load all libraries
invisible(lapply(required_packages, library, character.only = TRUE))


# -2- Ensure reproducibility with random sampling ---
set.seed(123)

# -3- Create folders ---
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



################################################################################################################
################################################################################################################



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
    # Call function to check for multicollinearity and perform stepwise BIC selection
    current_vars <- bioclim_selection(train_data)
    
  } else {
    # Collinearity check for manual variables
    message("Running collinearity check on bioclimatic variables...")
    cor_matrix <- cor(train_data[, current_vars], method = "spearman")
    
    print(round(cor_matrix, 2))
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
  
  
  # -1.3- Evaluate model ---
  
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
sp_predictors[[sp1]] <- c("bio3", "bio5", "bio8", "bio14")
sp_predictors[[sp2]] <- c( "bio6", "bio10", "bio12", "bio15")

# Task 1 Current SDM Results
sdm_results <- run_current_sdm(sp1, sp2, "Europe", sp_predictors)
current_maps <- sdm_results$maps
model_summary <- lapply(sdm_results$models, summary)

print(sdm_results$stats)
print(model_summary)



################################################################################################################
################################################################################################################



# ==============================================================================================================
#                                              --- TASK 2 ---
# ==============================================================================================================

# This section involves plotting the overlap in distribution of Species 1 and Species 2, and
# devising and calculating a metric for the degree of overlap between their ranges.
# Task 2 will be converted into one general function that can take any Species 1 and 2, and the
# results outputed from the previous task.

# The main steps of this function involves the following:
# -- 1: Store the binary map outputs from Task 1
# -- 2: Calculate the degree of overlap using principles from probability (Intersect/Union)
# -- 3: Plot the degree of overlap in region using colour-blind friendly palette.


# MAIN TASK 2 FUNCTION - DISTRIBUTION OVERLAP METRIC AND PLOT
# ==============================================================================================================
calculate_and_map_overlap <- function(species1, species2, sdm_results) {
  
  # -1- Calculating overlap of Species 1 and Species 2 using binary maps ---------------------------------------
  
  message("Calculating degree of overlap for ", species1, " and ", species2, "...")
  
  # -1.1- Extract binary rasters from Task 1 predict_and_map() function ---
  sp1_bin <- sdm_results$maps[[species1]]$binary_map
  sp2_bin <- sdm_results$maps[[species2]]$binary_map
  
  
  # -1.2- Calculate the degree of overlap (intersection/union) ---
  
  # Intersection (both species present)
  overlap_rast <- sp1_bin & sp2_bin
  
  # Union (at least one species is present)
  union_rast <- sp1_bin | sp2_bin
  
  # Count pixels using terra::global() - na.rm = TRUE ignores NAs/ocean points
  overlap_pixels <- terra::global(overlap_rast, "sum", na.rm = TRUE)[1, 1]
  union_pixels <- terra::global(union_rast, "sum", na.rm = TRUE)[1, 1]
  
  # Calculate percentage of overlap 
  overlap_metric <- (overlap_pixels / union_pixels) * 100
  
  cat("Degree of overlap of", species1, "and", species2, ":", round(overlap_metric, 2), "%\n")
  
  
  # -2- Plot the overlap ---------------------------------------------------------------------------------------
  
  # -2.1- Create a combined map for plotting ---
  # Binary raster only gives 0/1 - sp1 + (sp2 * 2) gives 4 distinct values
  
  # Possible values:
  # 0 = None (0 + 0)
  # 1 = Species 1 ONLY (1 + 0)
  # 2 = Species 2 ONLY (0 + 2)
  # 3 = BOTH (1 + 2)
  combined_map <- sp1_bin + (sp2_bin * 2)
  
  # Ensure areas with value 0 appear as no colour
  # combined_map[combined_map == 0] <- NA       <-- REMOVED 23/02/2026
  
  
  # -2.2- Plot the map ---
  # Define 4-colour okabe-ito palette assigned to values 0-3
  safe_palette <- c("grey80",   # 0
                    "#56B4E9",  # 1
                    "#D55E00",  # 2
                    "#009E73")  # 3
  
  plot(combined_map,
       main = paste("Distribution Overlap:", species1, "&", species2),
       col = safe_palette,
       legend = FALSE)
  
  # Add legend below the plot
  legend("bottom",
         inset = c(0, -0.20),
         legend = c("Neither",
                    paste(species1, "only"),
                    paste(species2, "only"),
                    "Overlap (Both)"),
         fill = safe_palette,
         bty = "n",
         cex = 0.8,
         xpd = TRUE,
         ncol = 2)
  
  message("Degree of overlap calculated and mapped successfully. Returning metrics...")
  return(list(
    metric = overlap_metric,
    overlap_raster = overlap_rast,
    combined_raster = combined_map))
  
}




# ==============================================================================================================
#                                               --- TASK 2 EXECUTION ---
# ==============================================================================================================

# Task 2 Distribution Overlap  Results
overlap_results <- calculate_and_map_overlap(sp1, sp2, sdm_results)



################################################################################################################
################################################################################################################



# ==============================================================================================================
#                                              --- TASK 3 ---
# ==============================================================================================================

# This section involves using a General Linear Model (GLM) to test whether the distribution
# of Species 1 at present is dependent on the distribution of Species 2.
# Task 3 will also be converted into a general function that takes the models and data
# output from Task 1 to execute the task for any 2 species.

# The main steps of this takes involves the following:
# -- 1: Take the training data for Species 1
# -- 2: Take continuous suitability map of Species 2 to extract probability of it
#       occurring at each Species 1 point.
# -- 3: Run glm() for Species 1 using bioclim variables AND Species 2 suitability.


# MAIN TASK 3 FUNCTION - DEPENDENCE OF SPECIES 1 ON SPECIES 2
# ==============================================================================================================
test_species_dependence <- function(species1, species2, sdm_results) {
  
  # -1- Extract Species 2 suitability values from Species 1 points ---------------------------------------------
  
  # Get dataset used for Species 1 model
  sp1_data <- sdm_results$data[[species1]]
  
  # Get cont. suitability map for Species 2
  sp2_map <- sdm_results$maps[[species2]]$continuous_map
  
  
  # -1.1- Extract suitability values ---
  # Create SpatVector for Species 1 points
  sp1_pts <- vect(sp1_data, geom = c("lon", "lat") , crs = crs(sp2_map))   # Ensure crs matches
  
  # Extract values
  sp2_vals <- terra::extract(sp2_map, sp1_pts)[, 2]
  
  # Add values as new predictor variable to Species 1's dataset
  sp1_data$sp2_suitability <- sp2_vals
  
  
  # -2- Fit the biological GLM ---------------------------------------------------------------------------------
  
  # Extract climate variables used in original Species 1 GLM
  sp1_vars <- attr(terms(sdm_results$models[[species1]]), "term.labels")
  
  # Dynamically build GLM formula: presence ~ biox + bioy + bioz + sp2_suitability
  formula_str <- paste("presence ~", paste(sp1_vars, collapse = " + "), "+ sp2_suitability")
  model_formula <- as.formula(formula_str)
  
  # Fit the model
  biotic_glm <- glm(model_formula, data = sp1_data, family = binomial)
  
  return(biotic_glm)
  
}




# ==============================================================================================================
#                                               --- TASK 3 EXECUTION ---
# ==============================================================================================================

# Task 3 Dependence Model Results
dependence_model <- test_species_dependence(sp1, sp2, sdm_results)

biotic_model_sum <- summary(dependence_model)
print(biotic_model_sum)



################################################################################################################
################################################################################################################



# ==============================================================================================================
#                                              --- TASK 4 ---
# ==============================================================================================================

# This section involves  predicting the future distribution of both Species 1 and Species 2
# separately using the CMIP6 data for future climate, and how the degree of overlap in ranges
# will change over time.
# Once again, this task will be functionalised so that predicting the future distribution
# for any two species can be done - it will utilise Task 1 (current SDM
# generation) outputs and the function for Task 2 (distribution overlap).

# The main steps of this function involves the following:
# -- 1: Download future climate data from CMIP6 and crop to user-chosen region.
# -- 2: Iterate future distribution prediction for each species.
# -- 3: Plot present vs. future suitability side-by-side.


# MAIN TASK 4 FUNCTION - FUTURE DISTRIBUTION PREDICTION
# ==============================================================================================================
run_future_sdm <- function(species1, species2, region, sdm_results) {
  
  # -0- Get region bounds ---
  if (!region %in% names(REGION_PRESETS)) {
    stop("Region not found. Please choose from the preset list.")
  }
  region_bounds <- REGION_PRESETS[[region]]
  
  
  # -1- Download and crop future climate data ------------------------------------------------------------------
  
  message("Downloading future climate data (CMIP6 CanESM5, SSP245, 2061-2080)...")
  
  future_clim <- cmip6_world(model = "CanESM5", var = "bio", ssp = "245", res = 10, 
                        time = "2061-2080", path = here("data/raw"))
  
  # Standardise names so future layers are bio1 - bio19 exactly
  names(future_clim) <- paste0("bio", 1:19)
  
  # Crop to region
  region_ext <- ext(region_bounds[1], region_bounds[2], region_bounds[3], region_bounds[4])
  future_crop <- crop(future_clim, region_ext)
  
  
  # -2- Iterate future distribution prediction for both species ------------------------------------------------
  
  future_maps <- list()
  
  for (sp_name in c(species1, species2)) {
    
    message(paste("Predicting future distribution for:", sp_name))
    
    # -2.1- Prepare future maps ---
    # Retrieve model, variables, and threshold from Task 1 model
    sp_model <- sdm_results$models[[sp_name]]
    sp_vars <- attr(terms(sp_model), "term.labels")
    sp_thr <- sdm_results$maps[[sp_name]]$threshold
    
    # Predict future cont. suitability
    future_pred <- terra::predict(future_crop[[sp_vars]], sp_model, type = "response")
    names(future_pred) <- "future_suitability"
    
    # Convert to future binary (presence/absence) map
    future_bin <- future_pred > sp_thr
    names(future_bin) <- "future_presence"
    
    
    # -2.2- Quantify current site range change ---
    # Create SpatVector for Species 1 point
    sp_data <- sdm_results$data[[sp_name]]
    sp_vect <- vect(sp_data, geom = c("lon", "lat"), crs = crs(future_pred))
    
    # Extract present and future values at current site
    cur_map <- sdm_results$maps[[sp_name]]$continuous_map
    cur_vals <- terra::extract(cur_map, sp_vect)[, 2]
    fu_vals <- terra::extract(future_pred, sp_vect)[, 2]
    
    cur_bin_site <- cur_vals >= sp_thr
    fu_bin_site <- fu_vals >= sp_thr
    
    range_loss <- sum(cur_bin_site & !fu_bin_site, na.rm = TRUE)
    range_gain <- sum(!cur_bin_site & fu_bin_site, na.rm = TRUE)
    
    cat("Current sites losing suitability:", range_loss, "\n")
    cat("Current sites gaining suitability:", range_gain, "\n\n")
    
    
    # -3- Plot comparison --------------------------------------------------------------------------------------
    
    par(mfrow = c(1, 2))
    plot(cur_map, main = paste("Present suitability:", sp_name))
    plot(future_pred, main = paste("Future suitability:", sp_name))
    
    # Store outputs
    future_maps[[sp_name]] <- list(
      continuous_map = future_pred,
      binary_map = future_bin)
    
    message("Completed generation of current and future SDMs for ", sp_name)
    cat("-----------------------------------------------------------------------\n")
    
  }
  
  # Reset plotting window
  par(mfrow = c(1, 1))
  
  return(list(maps = future_maps))
  
}




# ==============================================================================================================
#                                               --- TASK 4 EXECUTION ---
# ==============================================================================================================

# Run future SDMs
future_results <- run_future_sdm(sp1, sp2, "Europe", sdm_results)

# Re-calculate present overlap for comparison
overlap_present <- calculate_and_map_overlap(sp1, sp2, sdm_results)

# Calculate future overlap
overlap_future <- calculate_and_map_overlap(sp1, sp2, future_results)



################################################################################################################
################################################################################################################



# ==============================================================================================================
#                                              --- TASK 5 ---
# ==============================================================================================================

# This task will involve creating a pretty map function which can be used to map a
# publication-ready map to any of the suitability maps from the output.
# The user will be able to pick was CRS to use, what colour-blind friendly palette,
# and whether they want to add rivers/streams or not.


# MAIN TASK 5 FUNCTION - PUBLICATION-READY MAPS
# ==============================================================================================================
make_pretty_map <- function(suitability_map, 
                            species_name, 
                            target_crs, 
                            colour_palette, 
                            add_rivers,
                            file_name = "publication_map.pdf") {
  
  # -1- Prepare raster and country borders ---------------------------------------------------------------------
  
  map_proj <- project(suitability_map, target_crs)
  map_ext <- ext(map_proj)
  
  # Download contextual basemaps
  world_borders <- ne_countries(scale = 50, returnclass = "sf")
  world_borders_proj <- st_transform(world_borders, target_crs)
  
  
  # -2- Making pretty map --------------------------------------------------------------------------------------
  
  # Set up google font Montserrat for clean look
  font_add_google(name = "Montserrat", family = "Montserrat")
  showtext_auto()
  
  # Initialise ggplot with the raster
  pretty_map <- ggplot() +
    geom_sf(data = world_borders_proj, fill = NA, color = "grey90", linewidth = 0.3) +
    geom_spatraster(data = map_proj)
  
  # -2.1- Add rivers if TRUE ---
  if (add_rivers) {
    message("Downloading and adding river networks...")
    rivers <- ne_download(scale = 50, type = 'rivers_lake_centerlines',
                          category = 'physical', returnclass = "sf")
    rivers_proj <- st_transform(rivers, target_crs)
    pretty_map <- pretty_map + geom_sf(data = rivers_proj, color = "white",
                                       linewidth = 0.475, alpha = 0.7)
  }
  
  
  # -2.2- Add aesthetics ---
  pretty_map <- pretty_map +
    scale_fill_viridis_c(
      option = colour_palette,
      name = "Habitat\nSuitability",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      na.value = "transparent",
      guide = guide_colorbar(barwidth = 1.5, barheight = 15, title.position = "top")
    ) +
    # Add cartographic elements
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering()) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    # Use raster's region bounds
    coord_sf(xlim = c(map_ext[1], map_ext[2]),
             ylim = c(map_ext[3], map_ext[4]),
             expand = FALSE) +
    # -2.3- Add labels ---
    labs(
      title = bquote("Predicted distribution of" ~ italic(.(species_name))),
      subtitle = paste("Projection:", target_crs),
      x = "Longitude",
      y = "Latitude",
      caption = "Data acquired from GBIF and WorldClim"
    ) +
    theme_minimal(base_size = 12, base_family = "Montserrat") +
    theme(
      # Title/subtitle spacing
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40", margin = margin(b = 8)),
      plot.caption = element_text(size = 9, color = "grey50", hjust = 1, family = "Montserrat"),
      
      # Axis title spacing
      axis.title.x = element_text(margin = margin(t = 12), size = 11),
      axis.title.y = element_text(margin = margin(t = 12), size = 11),
      
      # Remove grid lines 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Legend styling
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold", family = "Montserrat"),
      legend.text = element_text(size = 9, family = "Montserrat"),
    )
  
  # -2.4- Save map automatically ---
  ggsave(here(file.path("outputs", "maps", file_name)), 
         pretty_map, 
         width = 7, height = 6, dpi = 300, bg = "white")
  
  return(pretty_map)
  
}




# ==============================================================================================================
#                                               --- TASK 5 EXECUTION ---
# ==============================================================================================================

# Store present cont. suitability map for your species 1
mistletoe_present_map <- sdm_results$maps[[sp1]]$continuous_map

# Run make_pretty_map() function
final_mistletoe_map <- make_pretty_map(
  mistletoe_present_map, 
  sp1,
  "EPSG:4326",
  "inferno",    # e.g. viridis, plasma, inferno
  add_rivers = TRUE,    # Set add_rivers to TRUE or FALSE for main river centrelines 
  "mistletoe_present_map.pdf")

# View map
print(final_mistletoe_map)


# Can do the exact same thing for species 2!
# Store present cont. suitability map for your species 2
oak_present_map <- sdm_results$maps[[sp2]]$continuous_map

final_oak_map <- make_pretty_map(
  oak_present_map, 
  sp2,
  "EPSG:4326",
  "inferno",
  add_rivers = TRUE,
  "oak_present_map.pdf")

print(final_oak_map)














