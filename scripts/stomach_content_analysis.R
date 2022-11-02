# =====================================================================
# Stomach Content Analysis
# Jonah Bacon
# 20 May 2022
# =====================================================================

# Load packages -----------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(ggsci)
library(vegan)

# Load data ---------------------------------------------------------------

stomach.df <- read.csv("data/Stomach_contents.csv", skip = 1, header = T)

stomach.df <- stomach.df %>% 
  separate(spp_ID, into = c("Species", "ID"), sep = "_", remove = FALSE)
stomach.df$Species <- as.factor(stomach.df$Species)
stomach.df$ID <- as.integer(stomach.df$ID)

prey.presence <- stomach.df %>% 
  select("spp_ID", "Species", "ID", "Prey_present", "Relative_fullness", "Total_contents_weight_g")

prey.count <- stomach.df %>% 
  select("spp_ID", "Species", "ID", contains(".Count"))
  
prey.weight <- stomach.df %>% 
  select("spp_ID", "Species", "ID", contains(".Weight"))

prey.percent <- stomach.df %>% 
  select("spp_ID", "Species", "ID", contains(".Relative_percent"))


# QAQC Prey data ----------------------------------------------------------

str(prey.presence)
unique(prey.presence$Species) # ARCS, LSCS, BDWF, HBWF
unique(prey.presence$Prey_present) # 0, 1
unique(prey.presence$Relative_fullness) # 0, 1, 2, 3
prey.presence$spp_ID[duplicated(prey.presence$spp_ID)]

str(prey.count)
unique(prey.count$Isopods.Count)
unique(prey.count$Amphipods.Count)
unique(prey.count$Bivalves.Count)
prey.count <- prey.count %>% 
  filter(Isopods.Count != "missing" & Amphipods.Count != "missing" & Bivalves.Count != "missing")
prey.count$Isopods.Count <- as.integer(prey.count$Isopods.Count)
prey.count$Amphipods.Count <- as.integer(prey.count$Amphipods.Count)
prey.count$Bivalves.Count <- as.integer(prey.count$Bivalves.Count)
str(prey.count)

str(prey.weight)
prey.weight$Vertebrate.Weight_g <- as.numeric(prey.weight$Vertebrate.Weight_g)
str(prey.weight)

str(prey.percent)

prey.count <- prey.count %>% 
  gather(key = "Prey_group", value = "Value", 4:16, na.rm = T) %>% 
  separate(Prey_group, into = c("Prey_group", "Measurement"))

prey.weight <- prey.weight %>% 
  gather(key = "Prey_group", value = "Value", 4:16, na.rm = T) %>% 
  separate(Prey_group, into = c("Prey_group", "Measurement"))
prey.weight$Measurement <- rep("Weight_g", length(prey.weight$Measurement))

prey.percent <- prey.percent %>% 
  gather(key = "Prey_group", value = "Value", 4:16, na.rm = T) %>% 
  separate(Prey_group, into = c("Prey_group", "Measurement"))
prey.percent$Measurement <- rep("Relative_percent", length(prey.percent$Measurement))


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


# Frequency of occurrence --------------------------------------------------

long.prey.df <- rbind(prey.count, prey.weight, prey.percent)

frequency.of.occurrence <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species) %>% 
  summarise(
    "Planerians" = sum(Prey_group == "Planerians"),
    "Isopods" = sum(Prey_group == "Isopods"),
    "Amphipods" = sum(Prey_group == "Amphipods"),
    "Bivalves" = sum(Prey_group == "Bivalves"),
    "Oligochaets" = sum(Prey_group == "Oligochaets"),
    "Insects" = sum(Prey_group == "Insects"),
    "Vertebrate" = sum(Prey_group == "Vertebrate"),
    "Chironomid" = sum(Prey_group == "Chironomid"),
    "Mysid" = sum(Prey_group == "Mysid"),
    "Fish" = sum(Prey_group == "Fish"),
    "Nematodes" = sum(Prey_group == "Nematodes"),
    "Unidentifiable" = sum(Prey_group == "Unidentifiable"),
    "Vegetation" = sum(Prey_group == "Vegetation")
  )
frequency.of.occurrence


# Rate of occurrence ------------------------------------------------------

rate.of.occurrence <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species) %>% 
  summarise(
    "Planerians" = round(sum(Prey_group == "Planerians")/length(unique(ID)),3),
    "Isopods" = round(sum(Prey_group == "Isopods")/length(unique(ID)),3),
    "Amphipods" = round(sum(Prey_group == "Amphipods")/length(unique(ID)),3),
    "Bivalves" = round(sum(Prey_group == "Bivalves")/length(unique(ID)),3),
    "Oligochaets" = round(sum(Prey_group == "Oligochaets")/length(unique(ID)),3),
    "Insects" = round(sum(Prey_group == "Insects")/length(unique(ID)),3),
    "Vertebrate" = round(sum(Prey_group == "Vertebrate")/length(unique(ID)),3),
    "Chironomid" = round(sum(Prey_group == "Chironomid")/length(unique(ID)),3),
    "Mysid" = round(sum(Prey_group == "Mysid")/length(unique(ID)),3),
    "Fish" = round(sum(Prey_group == "Fish")/length(unique(ID)),3),
    "Nematodes" = round(sum(Prey_group == "Nematodes")/length(unique(ID)),3),
    "Unidentifiable" = round(sum(Prey_group == "Unidentifiable")/length(unique(ID)),3),
    "Vegetation" = round(sum(Prey_group == "Vegetation")/length(unique(ID)),3)
  )
rate.of.occurrence


# Diversity indexes -------------------------------------------------------

diversity.index.counts <- long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(spp_ID, Species, ID) %>% 
  summarise(
    "Shannon" = diversity(Value),
    "Simpson" = diversity(Value, index = "simpson"),
    "InvSimpson" = diversity(Value, index = "invsimpson")
  )
diversity.index.counts

diversity.index.percents <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(spp_ID, Species, ID) %>% 
  summarise(
    "Shannon" = diversity(Value),
    "Simpson" = diversity(Value, index = "simpson"),
    "InvSimpson" = diversity(Value, index = "invsimpson")
  )
diversity.index.percents

diversity.index.weights <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(spp_ID, Species, ID) %>% 
  summarise(
    "Shannon" = diversity(Value),
    "Simpson" = diversity(Value, index = "simpson"),
    "InvSimpson" = diversity(Value, index = "invsimpson")
  )
diversity.index.weights

# Evenness index ----------------------------------------------------------

H.max.counts <- diversity.index.counts %>% 
  group_by(Species) %>% 
  summarise(
    "H.max" = max(Shannon)
  )
H.max.counts

evenness.index.counts <- diversity.index.counts %>% 
  summarise(
    "evenness" = ifelse(Species == "ARCS", Shannon/H.max.counts$H.max[1], 
                        ifelse(Species == "BDWF", Shannon/H.max.counts$H.max[2],
                               ifelse(Species == "HBWF", Shannon/H.max.counts$H.max[3],
                                      Shannon/H.max.counts$H.max[4]))) 
  )
evenness.index.counts

H.max.weights <- diversity.index.weights %>% 
  group_by(Species) %>% 
  summarise(
    "H.max" = max(Shannon)
  )
H.max.weights

evenness.index.weights <- diversity.index.weights %>% 
  summarise(
    "evenness" = ifelse(Species == "ARCS", Shannon/H.max.weights$H.max[1], 
                        ifelse(Species == "BDWF", Shannon/H.max.weights$H.max[2],
                               ifelse(Species == "HBWF", Shannon/H.max.weights$H.max[3],
                                      Shannon/H.max.weights$H.max[4]))) 
  )
evenness.index.weights

H.max.percents <- diversity.index.percents %>% 
  group_by(Species) %>% 
  summarise(
    "H.max" = max(Shannon)
  )
H.max.percents

evenness.index.percents <- diversity.index.percents %>% 
  summarise(
    "evenness" = ifelse(Species == "ARCS", Shannon/H.max.percents$H.max[1], 
                        ifelse(Species == "BDWF", Shannon/H.max.percents$H.max[2],
                               ifelse(Species == "HBWF", Shannon/H.max.percents$H.max[3],
                                      Shannon/H.max.percents$H.max[4]))) 
  )
evenness.index.percents


# Summarized diversity and evenness indices -------------------------------

summarized.diversity.evenness <- merge(diversity.index.counts, evenness.index.counts) %>% 
  group_by(Species) %>% 
  summarise(
    mean.Shannon = mean(Shannon),
    mean.Simpson = mean(Simpson),
    mean.invSimpson = mean(InvSimpson),
    mean.Evenness = mean(evenness),
  )
summarized.diversity.evenness

# Diet overlap ------------------------------------------------------------


# mean.prey.weight <- total.prey.weight %>% 
#   group_by(Species, Prey_group, Total_prey_weight) %>% 
#   summarise(
#     "N_stomachs_w_prey" = ifelse(Species == "ARCS", 46, 
#                                  ifelse(Species == "BDWF", 23,
#                                         ifelse(Species == "HBWF", 47,
#                                                56))),
#     "Mean_prey_weight" = Total_prey_weight/N_stomachs_w_prey
#   )
# mean.prey.weight



# total.prey.percent <- long.prey.df %>% 
#   filter(Measurement == "Relative_percent") %>% 
#   group_by(Species, Prey_group) %>% 
#   summarise(
#     "Total_prey_percent" = sum(Value)
#   )
# total.prey.percent
# 
# mean.prey.percent <- total.prey.percent %>% 
#   group_by(Species, Prey_group, Total_prey_percent) %>% 
#   summarise(
#     "N_stomachs_w_prey" = ifelse(Species == "ARCS", 46, 
#                                  ifelse(Species == "BDWF", 23,
#                                         ifelse(Species == "HBWF", 47,
#                                                56))),
#     "Mean_prey_percent" = Total_prey_percent/N_stomachs_w_prey
#   )
# mean.prey.percent


# Schoener index by prey count --------------------------------------------

long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(Species) %>% 
  summarise(
    sum(Value)
  )
# ARCS = 2451, BDWF = 1632, HBWF = 2842, LSCS = 4308

total.prey.count <- long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(Species, Prey_group) %>% 
  summarise(
    "Total_prey_count" = sum(Value),
    "Prey_count_percent" = unique(ifelse(Species == "ARCS", Total_prey_count/2451,
                                         ifelse(Species == "BDWF", Total_prey_count/1632,
                                                ifelse(Species == "HBWF", Total_prey_count/2842,
                                                       Total_prey_count/4308))))
  )
total.prey.count

ARCS.BDWF.count <- total.prey.count %>% 
  filter(Species == "ARCS" | Species == "BDWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.BDWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent[1] - Prey_percent[2]), Prey_percent)
  )
ARCS.HBWF.count <- total.prey.count %>% 
  filter(Species == "ARCS" | Species == "HBWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent[1] - Prey_percent[2]), Prey_percent)
  )
ARCS.LSCS.count <- total.prey.count %>% 
  filter(Species == "ARCS" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent[1] - Prey_percent[2]), Prey_percent)
  )
BDWF.HBWF.count <- total.prey.count %>% 
  filter(Species == "BDWF" | Species == "HBWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "BDWF.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent[1] - Prey_percent[2]), Prey_percent)
  )
BDWF.LSCS.count <- total.prey.count %>% 
  filter(Species == "BDWF" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "BDWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent[1] - Prey_percent[2]), Prey_percent)
  )
HBWF.LSCS.count <- total.prey.count %>% 
  filter(Species == "HBWF" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "HBWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent[1] - Prey_percent[2]), Prey_percent)
  )
prey.count.schoener <- merge(ARCS.BDWF, merge(ARCS.HBWF, merge(ARCS.LSCS, merge(BDWF.HBWF, merge(BDWF.LSCS, HBWF.LSCS, all = T), all = T), all = T), all = T), all = T)

schoener.index.count <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "schoener.index.count" = rep(NA,6))
schoener.index.count

i=1
for (i in 1:6) {
  schoener.index.count[i,2] <- 1 - 0.5*sum(prey.count.schoener[,i+1], na.rm = T)
}
schoener.index.count


# Schoener index by prey weight -------------------------------------------

test.df <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species) %>% 
  summarise(
    sum(Value)
  )
# ARCS = 16.7587, BDWF = 10.5100, HBWF = 54.5300, LSCS = 29.3190

total.prey.weight <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, Prey_group) %>% 
  summarise(
    "Total_prey_weight" = sum(Value),
    "Prey_weight_percent" = unique(ifelse(Species == "ARCS", Total_prey_weight/16.7587,
                                          ifelse(Species == "BDWF", Total_prey_weight/10.51,
                                                 ifelse(Species == "HBWF", Total_prey_weight/54.53,
                                                        Total_prey_weight/29.319))))
  )
total.prey.weight


ARCS.BDWF.weight <- total.prey.weight %>% 
  filter(Species == "ARCS" | Species == "BDWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.BDWF" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
ARCS.HBWF.weight <- total.prey.weight %>% 
  filter(Species == "ARCS" | Species == "HBWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
ARCS.LSCS.weight <- total.prey.weight %>% 
  filter(Species == "ARCS" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
BDWF.HBWF.weight <- total.prey.weight %>% 
  filter(Species == "BDWF" | Species == "HBWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "BDWF.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
BDWF.LSCS.weight <- total.prey.weight %>% 
  filter(Species == "BDWF" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "BDWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
HBWF.LSCS.weight <- total.prey.weight %>% 
  filter(Species == "HBWF" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "HBWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
prey.weight.schoener <- merge(ARCS.BDWF.weight, merge(ARCS.HBWF.weight, merge(ARCS.LSCS.weight, merge(BDWF.HBWF.weight, merge(BDWF.LSCS.weight, HBWF.LSCS.weight, all = T), all = T), all = T), all = T), all = T)

schoener.index.weight <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "schoener.index.weight" = rep(NA,6))
schoener.index.weight

i=1
for (i in 1:6) {
  schoener.index.weight[i,2] <- 1 - 0.5*sum(prey.weight.schoener[,i+1], na.rm = T)
}
schoener.index.weight


# Schoener index by prey relative percent ---------------------------------

long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(Species) %>% 
  summarise(
    sum(Value)
  )
# ARCS = 4501, BDWF = 2470, HBWF = 4577, LSCS = 5510

total.prey.percent <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(Species, Prey_group) %>% 
  summarise(
    "Total_prey_percent" = sum(Value),
    "Prey_percent_percent" = unique(ifelse(Species == "ARCS", Total_prey_percent/4501,
                                           ifelse(Species == "BDWF", Total_prey_percent/2470,
                                                  ifelse(Species == "HBWF", Total_prey_percent/4577,
                                                         Total_prey_percent/5510))))
  )
total.prey.percent

ARCS.BDWF.percent <- total.prey.percent %>% 
  filter(Species == "ARCS" | Species == "BDWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.BDWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
ARCS.HBWF.percent <- total.prey.percent %>% 
  filter(Species == "ARCS" | Species == "HBWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
ARCS.LSCS.percent <- total.prey.percent %>% 
  filter(Species == "ARCS" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
BDWF.HBWF.percent <- total.prey.percent %>% 
  filter(Species == "BDWF" | Species == "HBWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "BDWF.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
BDWF.LSCS.percent <- total.prey.percent %>% 
  filter(Species == "BDWF" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "BDWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
HBWF.LSCS.percent <- total.prey.percent %>% 
  filter(Species == "HBWF" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "HBWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
prey.percent.schoener <- merge(ARCS.BDWF.percent, merge(ARCS.HBWF.percent, merge(ARCS.LSCS.percent, merge(BDWF.HBWF.percent, merge(BDWF.LSCS.percent, HBWF.LSCS.percent, all = T), all = T), all = T), all = T), all = T)

schoener.index.percent <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "schoener.index.percent" = rep(NA,6))
schoener.index.percent

i=1
for (i in 1:6) {
  schoener.index.percent[i,2] <- 1 - 0.5*sum(prey.percent.schoener[,i+1], na.rm = T)
}
schoener.index.percent


# Combined Schoener index table -------------------------------------------

combined.schoener.index <- merge(schoener.index.count, merge(schoener.index.weight, schoener.index.percent))
combined.schoener.index


# ANOVA for differences in diet by predator -------------------------------

glm(predator ~ prey.count + prey.weight + prey.percent)
