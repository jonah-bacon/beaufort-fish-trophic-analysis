# =====================================================================
# Lens measurement analysis
# Jonah Bacon
# 18 October 2022
# =====================================================================


# Load libraries ----------------------------------------------------------

library(tidyverse)


# Load data ---------------------------------------------------------------

lens.meas <- read.csv("data/Lens_measurements.csv")
lgth.wt <- read.csv("data/Sample_LW.csv", colClasses = "character")

head(lens.meas)

length.df <- lgth.wt %>% 
  unite("Spp_ID", species:ID, sep= "_", remove = FALSE) %>% 
  select(Spp_ID, species, length_mm)

lens.dat <- lens.meas %>% 
  separate(Photo_ID, into = c("Spp_ID", "Photo_ID"), sep = 7) %>% 
  group_by(Spp_ID) %>% 
  summarise(diam_um = max(Measurement_um)) %>% 
  filter(diam_um > 0)

lens.length.df <- merge(length.df, lens.dat, by = "Spp_ID")
lens.length.df$length_mm <- as.integer(lens.length.df$length_mm)

ggplot(lens.length.df, aes(x = length_mm, y = diam_um, color = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, show.legend = T)

lm.coef.df <- lens.length.df %>% 
  group_by(species) %>% 
  summarise(
    slope = summary(lm(diam_um ~ length_mm))$coefficients[2],
    intercept = summary(lm(diam_um ~ length_mm))$coefficients[1],
    adj.r.sq = summary(lm(diam_um ~ length_mm))$adj.r.squared
  )

