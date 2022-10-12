# =====================================================================
# Stomach Content Analysis
# Jonah Bacon
# 20 May 2022
# =====================================================================

# Load packages -----------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(ggsci)


# Load data ---------------------------------------------------------------

stomach.df <- read.csv("data/stomach_contents_trial.csv")

prey.presence <- stomach.df %>% 
  select("spp_ID", "Prey_present", "Relative_fullness", "Total_contents_weight_g") %>% 
  filter(!row_number() %in% c(1)) %>% 
  separate(spp_ID, into = c("Species", "ID"), sep = "_", remove = FALSE)

unique(prey.presence$Species) # ARCS, LSCS, BDWF, HBWF
unique(prey.presence$Prey_present) # 0, 1
unique(prey.presence$Relative_fullness) # 0, 1, 2, 3

# Stomachs with prey ------------------------------------------------------

stomachs.w.prey <- prey.presence %>% 
  group_by(Species) %>% 
  summarise(
    "N_stomachs" = n(),
    "N_stomachs_w_prey" = sum(Prey_present > 0),
    "Proportion" = sum(Prey_present > 0)/length(Prey_present))
stomachs.w.prey


# Average fullness --------------------------------------------------------

average.fullness <- prey.presence %>% 
  filter(Prey_present == 1) %>% 
  group_by(Species) %>% 
  summarise("Stomachs<50%" = sum(Relative_fullness == 1),
            "Stomachs50-100%" = sum(Relative_fullness == 2),
            "Stomachs>100%" = sum(Relative_fullness == 3))
average.fullness


# Average stomach weight --------------------------------------------------

average.weight <- prey.presence %>% 
  filter(Prey_present == 1) %>% 
  group_by(Species) %>% 
  summarise("Average_weight_g" = mean(Total_contents_weight_g, na.rm = T))
average.weight


# Frequency of occurence --------------------------------------------------

planerians <- stomach.df %>% 
  select(spp_ID, Planerians.Flatworms, X, X.1) %>% 
  rename("Count" = "Planerians.Flatworms", "Weight_g" = "X", "Relative_percent" = X.1) %>% 
  filter(!row_number() %in% c(1)) %>% 
  mutate("Prey" = rep("planerians", length(spp_ID))) %>% 
  filter(Count > 0)
  
isopods <- stomach.df %>% 
  select(spp_ID, Isopods, X.2, X.3) %>% 
  rename("Count" = "Isopods", "Weight_g" = "X.2", "Relative_percent" = "X.3") %>% 
  filter(!row_number() %in% c(1)) %>% 
  mutate("Prey" = rep("isopods", length(spp_ID))) %>% 
  filter(Count > 0)

amphipods <- stomach.df %>% 
  select(spp_ID, Amphipods, X.4, X.5) %>% 
  rename("Count" = "Amphipods", "Weight_g" = "X.4", "Relative_percent" = "X.5") %>% 
  filter(!row_number() %in% c(1)) %>% 
  mutate("Prey" = rep("amphipods", length(spp_ID))) %>% 
  filter(Count > 0)

clams <- stomach.df %>% 
  select(spp_ID, Clams, X.6, X.7) %>% 
  rename("Count" = "Clams", "Weight_g" = "X.6", "Relative_percent" = "X.7") %>% 
  filter(!row_number() %in% c(1)) %>% 
  mutate("Prey" = rep("clams", length(spp_ID))) %>% 
  filter(Count > 0)

long.prey.df <- rbind(planerians, isopods, amphipods, clams)

frequency.of.occurrence <- long.prey.df %>% 
  separate(spp_ID, into = c("Species", "ID"), sep = "_", remove = FALSE) %>% 
  group_by(Species) %>% 
  summarise(
    "Planerians" = sum(Prey == "planerians")/35,
    "Isopods" = sum(Prey == "isopods")/35,
    "Amphipods" = sum(Prey == "amphipods")/35,
    "Clams" = sum(Prey == "clams")/35
  )
frequency.of.occurrence

