# Geo-data and Methods in R

## **Geospatial Species Distribution Modelling (SDM) Pipeline in R**
This repository contains a modular, reproducible pipeline for Species Distribution Modelling (SDM) in R. It was developed to analyse the spatial distribution, climatic drivers, and ecological interactions of two species (e.g., Loranthus europaeus and Quercus petraea), projecting their current and future habitat suitability under climate change scenarios.

The pipeline handles data acquisition (GBIF, WorldClim, CMIP6), data cleaning (spatial thinning, ocean filtering), statistical modeling (GLMs with stepwise BIC variable selection), and publication-ready cartography.

## **Installation**
To ensure reproducibility, this project uses the renv package to manage dependencies. **R Version: 4.5.2**

- _Option 1_: Restore Environment using renv (Recommended)
1. Clone the repository and open the R project.
2. R should automatically bootstrap renv.
3. If it doesn't, or to manually restore the exact package versions used in this project, run:
```r
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::restore()
```

- _Option 2:_ Manual Package Installation (Fallback)
If renv fails to restore the environment, you can use the built-in fallback script to install any missing dependencies globally:
```r
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
```

### **Quick Usage**
The main scripts are designed to be generalizable. You can run the entire pipeline for a new pair of interacting species simply by modifying the target variables at the top of the execution block, for example, in Task 1:
```r
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
```

### **Code structure:**
The codebase is abstracted into a series of helper and main execution functions to ensure modularity.


## **Core Data & Modeling Functions (Task 1):**
1. get_species_data(): Downloads occurrence limits from GBIF (rgbif), cleans coordinate data, removes spatial points located in the ocean, crops WorldClim data to the target region, and extracts climate values at valid occurrence points.
```r
get_species_data <- function(sp_name, region, region_bounds, ocean_vect, bioclim_global)
```

2. bioclim_selection(): An algorithmic variable selection tool. It removes highly correlated bioclimatic variables (Spearman > 0.7) and performs a Stepwise Generalized Linear Model (GLM) selection using the Bayesian Information Criterion (BIC) to find the most parsimonious predictors.
```r
bioclim_selection <- function(input_data)
```

3. fit_eval_glm(): Generates random background points (pseudo-absences), partitions data into training (70%) and testing (30%) sets, fits a binomial GLM, and evaluates model performance using Area Under the Curve (AUC).
```r
fit_eval_glm <- function(species_data, bioclim_crop, user_predictors = NULL)
```
_Note:_ If predictor_list = NULL is passed, the pipeline will automatically determine the most parsimonious bioclimatic variables using Spearman correlation thresholds and Stepwise BIC.

4. predict_and_map(): Projects the trained GLM onto the regional environmental rasters to generate continuous habitat suitability maps and calculates a biological threshold to generate binary Presence/Absence maps.
```r
predict_and_map <- function(sdm_model, bioclim_crop, current_vars, sp_name, species_data, eval_res)
```

5. run_current_sdm(): The main wrapper function that orchestrates the above four helpers for the target species pair, returning models, data frames, and mapped rasters.
```r
run_current_sdm <- function(species1, species2, region, predictor_list)
```

## **Analytical Functions (Tasks 2, 3 & 4):**
1. calculate_and_map_overlap(): Computes the geographic overlap between two species using an Intersection over Union (Jaccard) metric based on binary presence/absence maps. Outputs a 4-color combined discrete raster map.
```r
calculate_and_map_overlap <- function(species1, species2, sdm_results)
```

2. test_species_dependence(): Fits a "biotic" GLM to test if the presence of Species 1 is statistically dependent on the continuous habitat suitability of Species 2, while controlling for the original abiotic bioclimatic predictors.
```r
test_species_dependence <- function(species1, species2, sdm_results)
```

3. run_future_sdm(): Downloads CMIP6 future climate projections (e.g., CanESM5, SSP245), standardizes raster naming, predicts future habitat suitability, and calculates raw range loss and gain compared to present-day baselines.
```r
run_future_sdm <- function(species1, species2, region, sdm_results)
``` 

## **Publication ready map (Task 5):**
- make_pretty_map(): A comprehensive ggplot2 and tidyterra mapping tool. It standardizes projections, applies colorblind-friendly palettes (e.g., Viridis/Inferno), and adds cartographic elements such as scale bars, North arrows, and optional river networks for geographic context.
```r
make_pretty_map <- function(suitability_map, 
                            species_name, 
                            target_crs, 
                            colour_palette, 
                            add_rivers,
                            file_name = "publication_map.pdf")
```

