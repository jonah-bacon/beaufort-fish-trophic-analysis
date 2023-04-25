# =====================================================================
# Lens measurement analysis
# Jonah Bacon
# 27 February 2023
# =====================================================================


# Load libraries ----------------------------------------------------------

library(tidyverse)

# Load data ---------------------------------------------------------------

lens.df <- read.csv("data/data - Lens Measurements.csv")
lens.df <- lens.df %>% 
  filter(!Photo_ID %in% c("ARCS_133_03","ARCS_138_05")) %>% # Remove images that the photo malfunctioned
  mutate(
    Diameter_um = as.numeric(Measurement_um),
    Radius_um = Diameter_um/2)

lens.df %>%
  group_by(lens_ID) %>%
  filter(length(layer) != max(layer)) %>%
  summarise()
# 9 lenses had an image missed along the way (ARCS_133 and ARCS_138 had an image malfunction)

SIAposition <- lens.df %>% 
  filter(!is.na(Measurement_pixels)) %>% 
  group_by(lens_ID) %>% 
  summarize(
    Photo_ID = Photo_ID,
    SIAposition = ifelse(layer == max(layer), Radius_um/2, Radius_um[layer] - (Radius_um[layer]-Radius_um[layer + 1])/2),
    SIArange = ifelse(layer == max(layer), NA, Radius_um[layer]-Radius_um[layer + 1]))
SIAposition <- SIAposition[,-1]

lens.df <- lens.df %>% 
  left_join(SIAposition, by = "Photo_ID", keep = FALSE)

# Test for instances where the lens was measured larger after delamination occurred:
lens.df$test <- (lens.df$Radius_um - lens.df$SIAposition) < 0
lens.df %>% filter(test == TRUE) %>% summarise("N" = length(Photo_ID)) #31 instances where the lens was measured larger after peeling a lamina layer

save.image(file = "data/lens_measurements.RData")
# load("data/lens_measurements.RData")