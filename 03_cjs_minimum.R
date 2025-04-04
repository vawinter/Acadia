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
banding_data <- read_excel('MAPS atlantic forest BCR banding data.xlsx', sheet = 'MAPS_BANDING_capture_data')
station_data <- read.csv('MAPS atlantic forest BCR station info total.csv', header = T)

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
  distinct(STATION) %>%
  slice_sample(n = 6) %>%
  pull(STATION)

species_list <- c("BCCH", "BTNW", "HETH", "REVI")

# Then filter your data to only include those stations
full_data <- banding_data %>%
  filter(SPEC %in% species_list & STATION %in% sampled_stations) %>%
  select(STATION, year, SPEC, BAND)


#---------------------------------X
## Run this function for each species ----
#---------------------------------X
all_results <- list()
for (sp in species) {
  all_results[[sp]] <- test_min_recapture(full_data, sp)
}

