# =====================================================================
# BLE LTER baseline sample bulk SIA
# Jonah Bacon
# 10 April 2023
# =====================================================================

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggsci)
library(viridis)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

species.names <- c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")
names(species.names) <- c("ARCS","BDWF","HBWF","LSCS")


# Load data ---------------------------------------------------------------

biota <- read.csv("data/BLE_LTER_biota_stable_isotopes.csv", header = TRUE)

biota$node <- as.factor(biota$node)
biota$lagoon <- as.factor(biota$lagoon)
biota$location <- paste(biota$node, biota$lagoon, sep = ": ")
biota$station <- as.factor(biota$station)
biota$season <- as.factor(biota$season)
biota$taxa_id <- as.factor(biota$taxa_id)
biota$kingdom <- as.factor(biota$kingdom)
biota$phylum <- as.factor(biota$phylum)
biota$class <- as.factor(biota$class)
biota$order <- as.factor(biota$order)
biota$Family <- as.factor(biota$family)
biota$genus <- as.factor(biota$genus)
biota$species <- as.factor(biota$species)
biota$station_name <- as.factor(biota$station_name)
biota$habitat_type <- as.factor(biota$habitat_type)

# Summarize data ----------------------------------------------------------

df <- biota %>% 
  group_by(node, lagoon, Family) %>% 
  filter(!is.na(Family) & kingdom != "Chromista") %>% 
  summarise("mean.d13C" = mean(d13C, na.rm = T),
            "lower.d13C" = min(d13C, na.rm = T),
            "upper.d13C" = max(d13C, na.rm = T),
            "mean.d15N" = mean(d15N, na.rm = T),
            "lower.d15N" = min(d15N, na.rm = T),
            "upper.d15N" = max(d15N, na.rm = T))

biota %>% 
  group_by(node, Family) %>% 
  filter(!is.na(Family) & kingdom != "Chromista" & node == "Central") %>% 
  summarise("mean.d15N" = mean(d15N, na.rm = T))

# Plot data ---------------------------------------------------------------

ggplot(data = filter(df, node == "Central"), aes(x = mean.d13C, y = mean.d15N, color = Family)) +
  geom_point() +
  geom_linerange(aes(x = mean.d13C, ymin = lower.d15N, ymax = upper.d15N)) +
  geom_linerange(aes(y = mean.d15N, xmin = lower.d13C, xmax = upper.d13C)) +
  xlab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  facet_wrap(~lagoon, nrow = 1) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    # legend.position = "none",
    axis.line=element_line()
  )
# ggsave("figures/BLE.LTER.biota.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

