##########################################################################X
# CJS model to estimate survival from MAPS banding stations
#------------------------------------------------------------X
# What is the minimum number of birds per species necessary?
# Finding number by iterating through MAPS data and systematically
# removing samples each time.
# VAW
# 4/3/2025
##########################################################################X
# Clean env
rm(list=ls())
gc()
# reproducibility
set.seed(0235)

# Source dependencies
source("00_funs.R")

# Load banding data
banding_data <- read_excel('Data/MAPS atlantic forest BCR banding data.xlsx', sheet = 'MAPS_BANDING_capture_data')
station_data <- read.csv('Data/MAPS atlantic forest BCR station info total.csv', header = T)

############################################X
# Define NIMBLE Model ----
############################################X
cjs_code <- nimbleCode({
  for (s in 1:sp) {  # Loop over species
    for (l in 1:station) {  # Loop over stations
      # Global parameters for this species-station combo
      phi.mean[s, l] ~ dunif(0, 1)  # Mean survival probability
      p.mean[s, l] ~ dunif(0, 1)    # Mean capture probability
      
      for (t in 1:(time-1)) {
        phi[t, s, l] <- phi.mean[s, l]  # Use same survival for all time periods
      }
      
      for (t in 1:time) {
        p[t, s, l] <- p.mean[s, l]  # Use same capture prob for all time periods
      }
      
      for (i in 1:N) {  # Loop over individuals
        # First occasion - all individuals alive at first occasion (simplified)
        z[i, 1, s, l] <- 1
        y[i, 1, s, l] ~ dbern(p[1, s, l] * z[i, 1, s, l])
        
        for (t in 2:time) {
          # State process: survival from t-1 to t
          z[i, t, s, l] ~ dbern(phi[t-1, s, l] * z[i, t-1, s, l])
          
          # Observation process: detection at time t
          y[i, t, s, l] ~ dbern(p[t, s, l] * z[i, t, s, l])
        }
      }
    }
  }
})

############################################X
# Minimum captures? ----
############################################X
# First, get a random sample of 6 stations
sampled_stations <- banding_data %>%
  filter(SPEC %in% c("BCCH", "BTNW", "HETH", "REVI")) %>% 
  distinct(STATION) %>%
  slice_sample(n = 10) %>%
  pull(STATION)
#--------------------------------X
# Species set-up information ----
# c("BCCH", "BTNW", "HETH", "REVI")
#--------------------------------X
# BCCH (starting n: 64): iter = 350000, burn = 20000, t = 10; nmin = 18, 6 stations
# BTNW (starting n: 123):: iter = 350000, burn = 20000, t = 10; n min = 3 10 stations
# HETH (starting n: 90):: iter = 350000, burn = 20000, t = 10; n min = 2 6 stations
# REVI (starting n: 107):: iter = , burn = , t = ; n min = 3 6 stations

# SOI:
species <-  c("BCCH")
# filter data to only include those stations and species
full_data <- banding_data %>%
  filter(SPEC %in% species & STATION %in% sampled_stations) %>%
  select(STATION, year, SPEC, BAND)

#---------------------------------X
## Run this function for a single species ----
#---------------------------------X
test_min_recapture(full_data, species, iter = 350000, burn = 20000, t = 10)
# Check 6 sites,
# Check 1 site
#---------------------------------X
## Run this function for all species ----
#---------------------------------X
all_results <- lapply(species, function(sp) {
  cat("\nRunning for species:", sp, "\n")
  test_min_recapture(full_data, sp, iter = 350000, burn = 20000, t = 10)
})
names(all_results) <- species

# Create summary across all species
species_summary <- data.frame(
  Species = character(),
  Min_Sample_Size = numeric(),
  Mean_Precision = numeric(),
  Max_Rhat = numeric(),
  Min_ESS = numeric(),
  stringsAsFactors = FALSE
)

for (sp in names(results_list)) {
  if (!is.null(results_list[[sp]])) {
    min_n <- results_list[[sp]]$min_sample_size
    min_result <- results_list[[sp]]$all_results[[as.character(min_n)]]
    
    species_summary <- rbind(species_summary, data.frame(
      Species = sp,
      Min_Sample_Size = min_n,
      Mean_Precision = min_result$mean_precision,
      Max_Rhat = max(min_result$rhat, na.rm = TRUE),
      Min_ESS = min(min_result$ess, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
}

# Print final summary
cat("\n========================================\n")
cat("Final Summary Across Species\n")
cat("========================================\n")
print(species_summary)

# Save final results
saveRDS(full_data, file = "BCCH_full_data.rds")
saveRDS(samples, file = "BCCH_min_captures_2.rds")
write.csv(full_data, file = "species_min_sample_sizes.csv", row.names = FALSE)
