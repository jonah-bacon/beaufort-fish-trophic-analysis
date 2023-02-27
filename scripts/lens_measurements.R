# =====================================================================
# Lens measurement analysis
# Jonah Bacon
# 27 February 2023
# =====================================================================


# Load libraries ----------------------------------------------------------

library(tidyverse)

# Load data ---------------------------------------------------------------

lens.df <- read.csv("data/SIA_lens_measurements.csv")
lens.df$spp_ID <- substr(lens.df$Photo_ID, 1, 7)
lens.df$lens_ID <- sapply(strsplit(lens.df$Photo_ID, "_"), function(x) paste(x[1], x[2], sep = "_"))
lens.df$layer <- sapply(strsplit(lens.df$Photo_ID, "_"), function(x) paste(x[3]))
lens.df$layer <- as.integer(lens.df$layer)
lens.df$layer <- ifelse(is.na(lens.df$layer), 0, lens.df$layer)
lens.df$layer <- lens.df$layer+1
lens.df$Diameter_um <- as.numeric(lens.df$Measurement_um)
lens.df$Radius_um <- lens.df$Diameter_um/2

SIAposition <- lens.df %>% 
  group_by(lens_ID) %>% 
  summarise(
    Photo_ID = Photo_ID,
    SIAposition = ifelse(layer != max(layer), Radius_um[layer] - (Radius_um[layer]-Radius_um[layer + 1])/2, Radius_um[layer]/2))
SIAposition <- SIAposition[,-1]

lens.df <- lens.df %>% 
  left_join(SIAposition, by = "Photo_ID", keep = FALSE)

# Test for instances where the lens was measured larger after delamination occurred:
lens.df$test <- (lens.df$Radius_um - lens.df$SIAposition) < 0
lens.df %>% filter(test == TRUE) %>% summarise("N" = length(Photo_ID)) #37 instances

save.image(file = "data/lens_measurements.RData")
# load("data/lens_measurements.RData")

# lw.df<- read.csv("data/Sample_LW.csv", colClasses = "character")
# length.df <- lgth.wt %>% 
#   unite("Spp_ID", species:ID, sep= "_", remove = FALSE) %>% 
#   select(Spp_ID, species, length_mm)
# 
# lens.length.df <- merge(length.df, lens.dat, by = "Spp_ID")
# lens.length.df$length_mm <- as.integer(lens.length.df$length_mm)
# 
# ggplot(lens.length.df, aes(x = length_mm, y = diam_um, color = species)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F, show.legend = T)
# 
# lm.coef.df <- lens.length.df %>% 
#   group_by(species) %>% 
#   summarise(
#     slope = summary(lm(diam_um ~ length_mm))$coefficients[2],
#     intercept = summary(lm(diam_um ~ length_mm))$coefficients[1],
#     adj.r.sq = summary(lm(diam_um ~ length_mm))$adj.r.squared
#   )

