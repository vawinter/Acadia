##########################################################################X
# CJS model to estimate survival from MAPS banding stations
#------------------------------------------------------------X
# Fitting CJS model with MAPS
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
# Data prep ----
############################################X
# First, get a random sample of 6 stations
sampled_stations <- banding_data %>%
  distinct(STATION) %>%
  slice_sample(n = 6) %>%
  pull(STATION)

# Then filter your data to only include those stations
banding_data2 <- banding_data %>%
  filter(SPEC %in% species_list & STATION %in% sampled_stations) %>%
  select(STATION, year, SPEC, BAND)

# Filter for relevant species
species_list <- c('BCCH') #, 'BTNW', 'HETH', 'REVI')
filtered_data <- banding_data2 %>%
  filter(SPEC %in% species_list) %>%
  select(STATION, year, SPEC, BAND) 

# Assign numerical indices for categorical variables
filtered_data <- filtered_data %>%
  mutate(
    station_id = as.integer(as.factor(STATION)),
    species_id = as.integer(as.factor(SPEC)),
    year_id = as.integer(as.factor(year))
  )

# Get dataset dimensions
N <- length(unique(filtered_data$BAND))  # Total unique birds
sp <- length(species_list)  # Number of species
time <- length(unique(filtered_data$year))  # Number of years
station <- length(unique(filtered_data$STATION))  # Number of stations

# Create empty capture history array
capture_history <- array(0, dim = c(N, time, sp, station))

# Fill in capture history: 1 if the bird was captured, 0 otherwise
for (i in 1:nrow(filtered_data)) {
  bird <- which(unique(filtered_data$BAND) == filtered_data$BAND[i])
  yr <- filtered_data$year_id[i]
  stn <- filtered_data$station_id[i]
  spc <- filtered_data$species_id[i]
  
  capture_history[bird, yr, spc, stn] <- 1
}
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

# Define model constants
constants <- list(N = N, sp = sp, time = time, station = station)

# Prepare data list for NIMBLE
data_list <- list(y = capture_history)

# Initializing z based on capture history [N x time x station x sp]
z_inits <- init_z(capture_history)

# Define initial values
inits <- list(
  z = z_inits,
  phi.mean = array(0.5, dim = c(sp, station)),
  p.mean = array(0.5, dim = c(sp, station))
)

# Build and compile the model
cjs_model <- nimbleModel(
  cjs_code, 
  constants = constants, 
  data = data_list, 
  inits = inits
)

# Check model
cjs_model$calculate()  # Should not return -Inf if model is properly specified

# Compile model
compiled_model <- compileNimble(cjs_model)

# Set up MCMC
mcmc_conf <- configureMCMC(
  cjs_model, 
  monitors = c("phi", "p")
)

# Build and compile MCMC
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc)

# Run MCMC sampling
samples <- runMCMC(
  compiled_mcmc, 
  niter = 130000, 
  nburnin = 10000, 
  nchains = 2,
  thin = 10,
  samplesAsCodaMCMC = TRUE
)
#---------------------------------X
## Summarize results ----
#---------------------------------X
summary(samples)
MCMCvis::MCMCsummary(samples)
coda::raftery.diag(samples) # estimates how long the chain should be for reliable estimates

