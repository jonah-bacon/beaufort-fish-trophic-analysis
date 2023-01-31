# =====================================================================
# Stomach Content Analysis
# Jonah Bacon
# 20 May 2022
# =====================================================================

# Load packages -----------------------------------------------------------

library(ggsci)
library(vegan)
library(MASS)
library(ggrepel)
library(MicroNiche)
library(tidyverse)
library(glmmTMB)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

species <- c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")
names(species) <- c("ARCS","BDWF","HBWF","LSCS")

# Load data ---------------------------------------------------------------

stomach.df <- read.csv("data/Stomach_contents.csv", skip = 1, header = T)
# save(stomach.df, file = 
#        "data/stomach_content.Rdata")
# load(here::here("C:/Users/Bacon/Downloads/beaufort-fish-trophic-analysis/data/stomach_content.Rdata"))

stomach.df <- stomach.df %>% 
  separate(spp_ID, into = c("Species", "ID"), sep = "_", remove = FALSE)
stomach.df$Species <- as.factor(stomach.df$Species)
stomach.df$ID <- as.integer(stomach.df$ID)

lw.df <- read.csv("data/Sample_LW.csv", header = T)
lw.df <- lw.df %>% 
  rename(Species = species) %>% 
  unite(spp_ID, Species:ID, sep = "_", remove = FALSE)
lw.df$Species <- as.factor(lw.df$Species)
lw.df$date_collected <- as.Date(lw.df$date_collected, format = "%d-%b-%Y")
lw.df$year_collected <- format(lw.df$date_collected, "%Y")
lw.df$month_collected <- format(lw.df$date_collected, "%m")
lw.df$tperiod_collected <- 
  ifelse(between(lw.df$date_collected, as.Date("2021-07-01", format = "%Y-%m-%d"), as.Date("2021-07-15", format = "%Y-%m-%d")), 1,
  ifelse(between(lw.df$date_collected, as.Date("2022-07-01", format = "%Y-%m-%d"), as.Date("2022-07-15", format = "%Y-%m-%d")), 1,
  ifelse(between(lw.df$date_collected, as.Date("2021-07-16", format = "%Y-%m-%d"), as.Date("2021-07-31", format = "%Y-%m-%d")), 2,
  ifelse(between(lw.df$date_collected, as.Date("2022-07-16", format = "%Y-%m-%d"), as.Date("2022-07-31", format = "%Y-%m-%d")), 2,
  ifelse(between(lw.df$date_collected, as.Date("2021-08-01", format = "%Y-%m-%d"), as.Date("2021-08-15", format = "%Y-%m-%d")), 3,
  ifelse(between(lw.df$date_collected, as.Date("2022-08-01", format = "%Y-%m-%d"), as.Date("2022-08-15", format = "%Y-%m-%d")), 3, 4
  ))))))
lw.df$date_processed <- as.Date(lw.df$date_processed, format = "%d-%b-%Y")

sample.df <- lw.df %>% 
  select(spp_ID, date_collected, year_collected, month_collected, tperiod_collected, length_mm, weight_kg)

stomach.df <- merge(stomach.df, sample.df, by = "spp_ID")

prey.presence <- stomach.df %>% 
  select("spp_ID", "Species", "ID", "Prey_present", "Relative_fullness", "Total_contents_weight_g") %>% 
  merge(., sample.df, by = "spp_ID")

prey.count <- stomach.df %>% 
  select("spp_ID", "Species", "ID", contains(".Count")) %>% 
  merge(., sample.df, by = "spp_ID")

prey.weight <- stomach.df %>% 
  select("spp_ID", "Species", "ID", contains(".Weight")) %>% 
  merge(., sample.df, by = "spp_ID")

prey.percent <- stomach.df %>% 
  select("spp_ID", "Species", "ID", contains(".Relative_percent")) %>% 
  merge(., sample.df, by = "spp_ID")


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

prey.percent %>% 
  group_by(spp_ID) %>% 
  summarise(Total = sum(Value)) %>% 
  filter(Total != 100)

# Stomachs with prey ------------------------------------------------------

stomachs.w.prey <- prey.presence %>% 
  group_by(Species) %>% 
  summarise(
    "N_stomachs" = n(),
    "N_stomachs_w_prey" = sum(Prey_present > 0),
    "Proportion" = sum(Prey_present > 0)/length(Prey_present))
stomachs.w.prey

stomachs.w.prey %>% 
  gather(key = "Presence", value = "Value", 3:4, na.rm = T) %>% 
ggplot() +
  geom_col(aes(x = year_collected, y = Value, fill = Presence)) +
  facet_wrap(vars(Species)) +
  scale_fill_discrete(name = "Prey presence", labels = c("Absent", "Present"), guide = guide_legend(reverse = TRUE))

stomachs.w.prey %>% 
  gather(key = "Presence", value = "Value", 3:4, na.rm = T) %>% 
ggplot(aes(x = year_collected, y = Value, fill = Presence, label = Value)) +
  geom_col() +
  facet_wrap(vars(Species), labeller = labeller(Species = species)) +
  scale_fill_manual(values = cbPalette[c(1,2)], name = "Prey presence", labels = c("Absent", "Present"), guide = guide_legend(reverse = FALSE)) +
  ylab("Number of stomachs") +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  xlab("Year")

# Average fullness --------------------------------------------------------

average.fullness <- prey.presence %>% 
  filter(Prey_present == 1) %>% 
  group_by(Species, year_collected) %>% 
  summarise("Stomachs<50%" = sum(Relative_fullness == 1),
            "Stomachs50-100%" = sum(Relative_fullness == 2),
            "Stomachs>100%" = sum(Relative_fullness == 3))
average.fullness


# Average stomach weight --------------------------------------------------

average.weight <- prey.presence %>% 
  filter(Prey_present == 1) %>% 
  group_by(Species, year_collected) %>% 
  summarise("Average_weight_g" = mean(Total_contents_weight_g/weight_kg, na.rm = T))
average.weight

ggplot(prey.presence) +
  geom_boxplot(aes(x = year_collected, y = log(Total_contents_weight_g)), outlier.shape = NA) +
  geom_jitter(aes(x = year_collected, y = log(Total_contents_weight_g)), color = "black", alpha = 0.9) +
  facet_wrap(vars(Species))

ggplot(prey.presence) +
  geom_boxplot(aes(x = year_collected, y = Total_contents_weight_g), outlier.shape = NA) +
  geom_jitter(aes(x = year_collected, y = Total_contents_weight_g), color = "black", alpha = 0.9) +
  facet_wrap(vars(Species)) +
  scale_y_continuous(limits = c(0,5))

ggplot(prey.presence) +
  geom_histogram(aes(x = Total_contents_weight_g), col = "white") +
  facet_grid(rows = vars(Species), cols = vars(year_collected)) +
  scale_x_continuous(limits = c(0,5)) +
  scale_y_continuous(limits = c(0,14), breaks = seq(0,14,2))

# Frequency of occurrence --------------------------------------------------

long.prey.df <- rbind(prey.count, prey.weight, prey.percent)
long.prey.df$lg_group = floor(long.prey.df$length_mm/50)

frequency.of.occurrence <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected) %>% 
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
    "n.stomachs" = length(unique(ID)),
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


# Comparison between years ------------------------------------------------

## Weight
df1 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected, lg_group) %>% 
  summarise(
    "Total.Stomach.Wt" = sum(Value)
    )
df2 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected, Prey_group, lg_group) %>% 
  summarise(
    "Total.Prey.Wt" = sum(Value)
  )
df3 <- merge(df1,df2, by = c("Species", "year_collected", "lg_group"))
df3$Wt.percent <- df3$Total.Prey.Wt/df3$Total.Stomach.Wt

df13 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    "Total.Stomach.Wt" = sum(Value)
  )
df14 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected, Prey_group) %>% 
  summarise(
    "Total.Prey.Wt" = sum(Value)
  )
df15 <- merge(df13,df14, by = c("Species", "year_collected"))
df15$Wt.percent <- df15$Total.Prey.Wt/df15$Total.Stomach.Wt

ggplot(df3, aes(x = year_collected, y = Wt.percent)) +
  geom_col(aes(fill = Prey_group)) +
  facet_wrap(vars(Species), nrow = 2)


## Count
df4 <- long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(Species, year_collected, lg_group) %>% 
  summarise(
    "Total.Stomach.Count" = sum(Value)
  )
df5 <- long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(Species, year_collected, Prey_group, lg_group) %>% 
  summarise(
    "Total.Prey.Count" = sum(Value)
  )
df6 <- merge(df4,df5, by = c("Species", "year_collected", "lg_group"))
df6$Count.percent <- df6$Total.Prey.Count/df6$Total.Stomach.Count

df16 <- long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    "Total.Stomach.Count" = sum(Value)
  )
df17 <- long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(Species, year_collected, Prey_group) %>% 
  summarise(
    "Total.Prey.Count" = sum(Value)
  )
df18 <- merge(df16,df17, by = c("Species", "year_collected"))
df18$Count.percent <- df18$Total.Prey.Count/df18$Total.Stomach.Count

ggplot(df6, aes(x = year_collected, y = Count.percent)) +
  geom_col(aes(fill = Prey_group)) +
  facet_wrap(vars(Species), nrow = 2)


## Occurrence
df7 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected, lg_group) %>% 
  summarise(
    "Total_N_stomachs_w_prey" = length(unique(spp_ID))
  )
df8 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected, Prey_group, lg_group) %>% 
  summarise(
    "N_stomachs_w_prey" = length(unique(spp_ID))
  )
df9 <- merge(df7,df8, by = c("Species", "year_collected", "lg_group"))
df9$Occurrence.percent <- df9$N_stomachs_w_prey/df9$Total_N_stomachs_w_prey

df19 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    "Total_N_stomachs_w_prey" = length(unique(spp_ID))
  )
df20 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected, Prey_group) %>% 
  summarise(
    "N_stomachs_w_prey" = length(unique(spp_ID))
  )
df21 <- merge(df19,df20, by = c("Species", "year_collected"))
df21$Occurrence.percent <- df21$N_stomachs_w_prey/df21$Total_N_stomachs_w_prey

ggplot(df9, aes(x = year_collected, y = Occurrence.percent)) +
  geom_col(aes(fill = Prey_group)) +
  facet_wrap(vars(Species), nrow = 2)


## Relative percent by volume
df10 <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(Species, year_collected, lg_group) %>% 
  summarise(
    "Total.Stomach.Percent" = sum(Value)
  )
df11 <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(Species, year_collected, Prey_group, lg_group) %>% 
  summarise(
    "Total.Prey.Percent" = sum(Value)
  )
df12 <- merge(df10,df11, by = c("Species", "year_collected", "lg_group"))
df12$Percent.percent <- df12$Total.Prey.Percent/df12$Total.Stomach.Percent
ggplot(df12, aes(x = year_collected, y = Percent.percent)) +
  geom_col(aes(fill = Prey_group)) +
  facet_wrap(vars(Species), nrow = 2)

df22 <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    "Total.Stomach.Percent" = sum(Value)
  )
df23 <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(Species, year_collected, Prey_group
           ) %>% 
  summarise(
    "Total.Prey.Percent" = sum(Value)
  )
df24 <- merge(df22,df23, by = c("Species", "year_collected"))
df24$Percent.percent <- df24$Total.Prey.Percent/df24$Total.Stomach.Percent

ggplot(df12, aes(x = year_collected, y = Percent.percent)) +
  geom_col(aes(fill = Prey_group)) +
  facet_wrap(vars(Species), nrow = 2)



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

diversity.index.counts.year <- long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(spp_ID, Species, ID, year_collected) %>% 
  summarise(
    "Shannon" = diversity(Value),
    "Simpson" = diversity(Value, index = "simpson"),
    "InvSimpson" = diversity(Value, index = "invsimpson")
  )
diversity.index.counts.year

diversity.index.percents <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(spp_ID, Species, ID) %>% 
  summarise(
    "Shannon" = diversity(Value),
    "Simpson" = diversity(Value, index = "simpson"),
    "InvSimpson" = diversity(Value, index = "invsimpson")
  )
diversity.index.percents

diversity.index.percents.year <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(spp_ID, Species, ID, year_collected) %>% 
  summarise(
    "Shannon" = diversity(Value),
    "Simpson" = diversity(Value, index = "simpson"),
    "InvSimpson" = diversity(Value, index = "invsimpson")
  )
diversity.index.percents.year

diversity.index.weights <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(spp_ID, Species, ID) %>% 
  summarise(
    "Shannon" = diversity(Value),
    "Simpson" = diversity(Value, index = "simpson"),
    "InvSimpson" = diversity(Value, index = "invsimpson")
  )
diversity.index.weights

diversity.index.weights.year <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(spp_ID, Species, ID, year_collected) %>% 
  summarise(
    "Shannon" = diversity(Value),
    "Simpson" = diversity(Value, index = "simpson"),
    "InvSimpson" = diversity(Value, index = "invsimpson")
  )
diversity.index.weights.year

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


## By year:

H.max.counts.year <- diversity.index.counts.year %>% 
  group_by(Species,year_collected) %>% 
  summarise(
    "H.max" = max(Shannon)
  )
H.max.counts.year

evenness.index.counts.year <- diversity.index.counts.year %>% 
  summarise(
    "evenness" = ifelse(Species == "ARCS" & year_collected == "2021", Shannon/H.max.counts.year$H.max[1], 
                  ifelse(Species == "ARCS" & year_collected == "2022", Shannon/H.max.counts.year$H.max[2],
                    ifelse(Species == "BDWF"  & year_collected == "2021", Shannon/H.max.counts.year$H.max[3],
                      ifelse(Species == "BDWF"  & year_collected == "2022", Shannon/H.max.counts.year$H.max[4],
                         ifelse(Species == "HBWF"  & year_collected == "2021", Shannon/H.max.counts.year$H.max[5],
                           ifelse(Species == "HBWF"  & year_collected == "2022", Shannon/H.max.counts.year$H.max[6],
                                      Shannon/H.max.counts.year$H.max[7])))))) 
  )
evenness.index.counts.year

H.max.weights.year <- diversity.index.weights.year %>% 
  group_by(Species,year_collected) %>% 
  summarise(
    "H.max" = max(Shannon)
  )
H.max.weights.year

evenness.index.weights.year <- diversity.index.weights.year %>% 
  summarise(
    "evenness" = ifelse(Species == "ARCS" & year_collected == "2021", Shannon/H.max.weights.year$H.max[1], 
                        ifelse(Species == "ARCS" & year_collected == "2022", Shannon/H.max.weights.year$H.max[2],
                               ifelse(Species == "BDWF"  & year_collected == "2021", Shannon/H.max.weights.year$H.max[3],
                                      ifelse(Species == "BDWF"  & year_collected == "2022", Shannon/H.max.weights.year$H.max[4],
                                             ifelse(Species == "HBWF"  & year_collected == "2021", Shannon/H.max.weights.year$H.max[5],
                                                    ifelse(Species == "HBWF"  & year_collected == "2022", Shannon/H.max.weights.year$H.max[6],
                                                           Shannon/H.max.weights.year$H.max[7])))))) 
  )
evenness.index.weights.year

H.max.percents.year <- diversity.index.percents.year %>% 
  group_by(Species,year_collected) %>% 
  summarise(
    "H.max" = max(Shannon)
  )
H.max.percents.year

evenness.index.percents.year <- diversity.index.percents.year %>% 
  summarise(
    "evenness" = ifelse(Species == "ARCS" & year_collected == "2021", Shannon/H.max.percents.year$H.max[1], 
                        ifelse(Species == "ARCS" & year_collected == "2022", Shannon/H.max.percents.year$H.max[2],
                               ifelse(Species == "BDWF"  & year_collected == "2021", Shannon/H.max.percents.year$H.max[3],
                                      ifelse(Species == "BDWF"  & year_collected == "2022", Shannon/H.max.percents.year$H.max[4],
                                             ifelse(Species == "HBWF"  & year_collected == "2021", Shannon/H.max.percents.year$H.max[5],
                                                    ifelse(Species == "HBWF"  & year_collected == "2022", Shannon/H.max.percents.year$H.max[6],
                                                           Shannon/H.max.percents.year$H.max[7])))))) 
  )
evenness.index.percents.year

# Summarized diversity and evenness indices -------------------------------

summarized.diversity.evenness.counts <- merge(diversity.index.counts, evenness.index.counts) %>% 
  group_by(Species) %>% 
  summarise(
    mean.Shannon = mean(Shannon),
    mean.Simpson = mean(Simpson),
    mean.invSimpson = mean(InvSimpson),
    mean.Evenness = mean(evenness),
  )
summarized.diversity.evenness.counts

summarized.diversity.evenness.weights <- merge(diversity.index.weights, evenness.index.weights) %>% 
  group_by(Species) %>% 
  summarise(
    mean.Shannon = mean(Shannon),
    mean.Simpson = mean(Simpson),
    mean.invSimpson = mean(InvSimpson),
    mean.Evenness = mean(evenness),
  )
summarized.diversity.evenness.weights

summarized.diversity.evenness.percents <- merge(diversity.index.percents, evenness.index.percents) %>% 
  group_by(Species) %>% 
  summarise(
    mean.Shannon = mean(Shannon),
    mean.Simpson = mean(Simpson),
    mean.invSimpson = mean(InvSimpson),
    mean.Evenness = mean(evenness),
  )
summarized.diversity.evenness.percents

summarized.diversity.evenness.counts.year <- merge(diversity.index.counts.year, evenness.index.counts.year) %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    mean.Shannon = mean(Shannon),
    mean.Simpson = mean(Simpson),
    mean.invSimpson = mean(InvSimpson),
    mean.Evenness = mean(evenness),
  )
summarized.diversity.evenness.counts.year

summarized.diversity.evenness.weights.year <- merge(diversity.index.weights.year, evenness.index.weights.year) %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    mean.Shannon = mean(Shannon),
    mean.Simpson = mean(Simpson),
    mean.invSimpson = mean(InvSimpson),
    mean.Evenness = mean(evenness),
  )
summarized.diversity.evenness.weights.year

summarized.diversity.evenness.percents.year <- merge(diversity.index.percents.year, evenness.index.percents.year) %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    mean.Shannon = mean(Shannon),
    mean.Simpson = mean(Simpson),
    mean.invSimpson = mean(InvSimpson),
    mean.Evenness = mean(evenness),
  )
summarized.diversity.evenness.percents.year

# Diet overlap ------------------------------------------------------------

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
    "ARCS.BDWF" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
ARCS.HBWF.count <- total.prey.count %>% 
  filter(Species == "ARCS" | Species == "HBWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
ARCS.LSCS.count <- total.prey.count %>% 
  filter(Species == "ARCS" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "ARCS.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
BDWF.HBWF.count <- total.prey.count %>% 
  filter(Species == "BDWF" | Species == "HBWF") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "BDWF.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
BDWF.LSCS.count <- total.prey.count %>% 
  filter(Species == "BDWF" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "BDWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
HBWF.LSCS.count <- total.prey.count %>% 
  filter(Species == "HBWF" | Species == "LSCS") %>% 
  group_by(Prey_group) %>% 
  summarise(
    "HBWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
prey.count.schoener <- merge(ARCS.BDWF.count, merge(ARCS.HBWF.count, merge(ARCS.LSCS.count, merge(BDWF.HBWF.count, merge(BDWF.LSCS.count, HBWF.LSCS.count, all = T), all = T), all = T), all = T), all = T)

schoener.index.count <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "schoener.index.count" = rep(NA,6))
schoener.index.count

i=1
for (i in 1:6) {
  schoener.index.count[i,2] <- 1 - 0.5*sum(prey.count.schoener[,i+1], na.rm = T)
}
schoener.index.count


# Schoener index by prey weight -------------------------------------------

long.prey.df %>% 
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




# Schoener index by prey count by year --------------------------------------------

long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    sum(Value)
  )
# ARCS 2021 = 2, ARCS 2022 = 2449, BDWF 2021 = 3, BDWF 2022 = 1614, HBWF 2021 = 225, HBWF 2022 = 2617, LSCS 2022 = 4308

total.prey.count2 <- long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  group_by(Species, Prey_group, year_collected) %>% 
  summarise(
    "Total_prey_count" = sum(Value),
    "Prey_count_percent" = unique(ifelse(Species == "ARCS" & year_collected == 2021, Total_prey_count/2,
                                  ifelse(Species == "ARCS" & year_collected == 2022, Total_prey_count/2449,
                                  ifelse(Species == "BDWF" & year_collected == 2021, Total_prey_count/3,
                                  ifelse(Species == "BDWF" & year_collected == 2022, Total_prey_count/1614,
                                  ifelse(Species == "HBWF" & year_collected == 2021, Total_prey_count/225,
                                  ifelse(Species == "HBWF" & year_collected == 2022, Total_prey_count/2617,
                                     Total_prey_count/4308)))))))
  )
total.prey.count2

ARCS.BDWF.count2 <- total.prey.count2 %>% 
  filter(Species == "ARCS" | Species == "BDWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.BDWF" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
ARCS.HBWF.count2 <- total.prey.count2 %>% 
  filter(Species == "ARCS" | Species == "HBWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
ARCS.LSCS.count2 <- total.prey.count2 %>% 
  filter(Species == "ARCS" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
BDWF.HBWF.count2 <- total.prey.count2 %>% 
  filter(Species == "BDWF" | Species == "HBWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "BDWF.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
BDWF.LSCS.count2 <- total.prey.count2 %>% 
  filter(Species == "BDWF" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "BDWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
HBWF.LSCS.count2 <- total.prey.count2 %>% 
  filter(Species == "HBWF" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "HBWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_count_percent[1] - Prey_count_percent[2]), Prey_count_percent)
  )
prey.count.schoener2 <- merge(ARCS.BDWF.count2, merge(ARCS.HBWF.count2, merge(ARCS.LSCS.count2, merge(BDWF.HBWF.count2, merge(BDWF.LSCS.count2, HBWF.LSCS.count2, all = T), all = T), all = T), all = T), all = T)

schoener.index.count2021 <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "year" = c(rep("2021", 6)), "schoener.index.count" = rep(NA,6))
schoener.index.count2021
schoener.index.count2022 <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "year" = c(rep("2022", 6)), "schoener.index.count" = rep(NA,6))
schoener.index.count2022

i=1
for (i in 1:6) {
  schoener.index.count2021[i,3] <- 1 - 0.5*sum(filter(prey.count.schoener2, year_collected == "2021")[,i+2], na.rm = T)
}
schoener.index.count2021
schoener.index.count2021 <- schoener.index.count2021[-c(3,5,6),]

i=1
for (i in 1:6) {
  schoener.index.count2022[i,3] <- 1 - 0.5*sum(filter(prey.count.schoener2, year_collected == "2022")[,i+2], na.rm = T)
}
schoener.index.count2022
schoener.index.count.year <- rbind(schoener.index.count2021, schoener.index.count2022)

# Schoener index by prey weight by year -------------------------------------------

long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    round(sum(Value), 4)
  )
# 

total.prey.weight2 <- long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  group_by(Species, Prey_group, year_collected) %>% 
  summarise(
    "Total_prey_weight" = sum(Value),
    "Prey_weight_percent" = unique(ifelse(Species == "ARCS" & year_collected == 2021, Total_prey_weight/0.1460,
                                         ifelse(Species == "ARCS" & year_collected == 2022, Total_prey_weight/16.613,
                                                ifelse(Species == "BDWF" & year_collected == 2021, Total_prey_weight/0.2030,
                                                       ifelse(Species == "BDWF" & year_collected == 2022, Total_prey_weight/10.526,
                                                              ifelse(Species == "HBWF" & year_collected == 2021, Total_prey_weight/12.958,
                                                                     ifelse(Species == "HBWF" & year_collected == 2022, Total_prey_weight/49.795,
                                                                            Total_prey_weight/30.013)))))))
  )
total.prey.weight2

ARCS.BDWF.weight2 <- total.prey.weight2 %>% 
  filter(Species == "ARCS" | Species == "BDWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.BDWF" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
ARCS.HBWF.weight2 <- total.prey.weight2 %>% 
  filter(Species == "ARCS" | Species == "HBWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
ARCS.LSCS.weight2 <- total.prey.weight2 %>% 
  filter(Species == "ARCS" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
BDWF.HBWF.weight2 <- total.prey.weight2 %>% 
  filter(Species == "BDWF" | Species == "HBWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "BDWF.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
BDWF.LSCS.weight2 <- total.prey.weight2 %>% 
  filter(Species == "BDWF" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "BDWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
HBWF.LSCS.weight2 <- total.prey.weight2 %>% 
  filter(Species == "HBWF" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "HBWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_weight_percent[1] - Prey_weight_percent[2]), Prey_weight_percent)
  )
prey.weight.schoener2 <- merge(ARCS.BDWF.weight2, merge(ARCS.HBWF.weight2, merge(ARCS.LSCS.weight2, merge(BDWF.HBWF.weight2, merge(BDWF.LSCS.weight2, HBWF.LSCS.weight2, all = T), all = T), all = T), all = T), all = T)

schoener.index.weight2021 <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "year" = c(rep("2021", 6)), "schoener.index.weight" = rep(NA,6))
schoener.index.weight2021
schoener.index.weight2022 <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "year" = c(rep("2022", 6)), "schoener.index.weight" = rep(NA,6))
schoener.index.weight2022

i=1
for (i in 1:6) {
  schoener.index.weight2021[i,3] <- 1 - 0.5*sum(filter(prey.weight.schoener2, year_collected == "2021")[,i+2], na.rm = T)
}
schoener.index.weight2021
schoener.index.weight2021 <- schoener.index.weight2021[-c(3,5,6),]

i=1
for (i in 1:6) {
  schoener.index.weight2022[i,3] <- 1 - 0.5*sum(filter(prey.weight.schoener2, year_collected == "2022")[,i+2], na.rm = T)
}
schoener.index.weight2022
schoener.index.weight.year <- rbind(schoener.index.weight2021, schoener.index.weight2022)



# Schoener index by prey relative percent by year ---------------------------------

long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(Species, year_collected) %>% 
  summarise(
    round(sum(Value), 4)
  )
# 

total.prey.percent2 <- long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  group_by(Species, Prey_group, year_collected) %>% 
  summarise(
    "Total_prey_percent" = sum(Value),
    "Prey_percent_percent" = unique(ifelse(Species == "ARCS" & year_collected == 2021, Total_prey_percent/100,
                                          ifelse(Species == "ARCS" & year_collected == 2022, Total_prey_percent/4400,
                                                 ifelse(Species == "BDWF" & year_collected == 2021, Total_prey_percent/200,
                                                        ifelse(Species == "BDWF" & year_collected == 2022, Total_prey_percent/2300,
                                                               ifelse(Species == "HBWF" & year_collected == 2021, Total_prey_percent/1500,
                                                                      ifelse(Species == "HBWF" & year_collected == 2022, Total_prey_percent/3100,
                                                                             Total_prey_percent/5600)))))))
  )
total.prey.percent2

ARCS.BDWF.percent2 <- total.prey.percent2 %>% 
  filter(Species == "ARCS" | Species == "BDWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.BDWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
ARCS.HBWF.percent2 <- total.prey.percent2 %>% 
  filter(Species == "ARCS" | Species == "HBWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
ARCS.LSCS.percent2 <- total.prey.percent2 %>% 
  filter(Species == "ARCS" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "ARCS.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
BDWF.HBWF.percent2 <- total.prey.percent2 %>% 
  filter(Species == "BDWF" | Species == "HBWF") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "BDWF.HBWF" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
BDWF.LSCS.percent2 <- total.prey.percent2 %>% 
  filter(Species == "BDWF" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "BDWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
HBWF.LSCS.percent2 <- total.prey.percent2 %>% 
  filter(Species == "HBWF" | Species == "LSCS") %>% 
  group_by(Prey_group, year_collected) %>% 
  summarise(
    "HBWF.LSCS" = ifelse(length(Prey_group) == 2, abs(Prey_percent_percent[1] - Prey_percent_percent[2]), Prey_percent_percent)
  )
prey.percent.schoener2 <- merge(ARCS.BDWF.percent2, merge(ARCS.HBWF.percent2, merge(ARCS.LSCS.percent2, merge(BDWF.HBWF.percent2, merge(BDWF.LSCS.percent2, HBWF.LSCS.percent2, all = T), all = T), all = T), all = T), all = T)

schoener.index.percent2021 <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "year" = c(rep("2021", 6)), "schoener.index.percent" = rep(NA,6))
schoener.index.percent2021
schoener.index.percent2022 <- data.frame("Species.interactions" = c("ARCS:BDWF", "ARCS:HBWF","ARCS:LSCS","BDWF:HBWF","BDWF:LSCS", "HBWF:LSCS"), "year" = c(rep("2022", 6)), "schoener.index.percent" = rep(NA,6))
schoener.index.percent2022

i=1
for (i in 1:6) {
  schoener.index.percent2021[i,3] <- 1 - 0.5*sum(filter(prey.percent.schoener2, year_collected == "2021")[,i+2], na.rm = T)
}
schoener.index.percent2021
schoener.index.percent2021 <- schoener.index.percent2021[-c(3,5,6),]

i=1
for (i in 1:6) {
  schoener.index.percent2022[i,3] <- 1 - 0.5*sum(filter(prey.percent.schoener2, year_collected == "2022")[,i+2], na.rm = T)
}
schoener.index.percent2022
schoener.index.percent.year <- rbind(schoener.index.percent2021, schoener.index.percent2022)


# Combined Schoener index table by year -------------------------------------------

combined.schoener.index.year <- merge(schoener.index.count.year, merge(schoener.index.weight.year, schoener.index.percent.year))
combined.schoener.index.year




# Levin's niche breadth ---------------------------------------------------

count.dat <- filter(long.prey.df, Measurement == "Count" & spp_ID != "ARCS_106")
levins.df <- data.frame(count.dat %>% 
  pivot_wider(id_cols = spp_ID, names_from = Prey_group, values_from = Value))
# levins.df <- data.frame(count.dat %>% 
#                           pivot_longer(id_cols = Species, names_from = Prey_group, values_from = sum(Value)))
levins.df[is.na(levins.df)] <- 0
levins.sampleInfo <- as.character(count.dat$Species)

levins.Bn(levins.df, 11, levins.sampleInfo)

# Save/load RData file ---------------------------------------------------------

save.image(file = "data/stomach_content_workspace.RData")
load("data/stomach_content_workspace.RData")

# Boxplots ----------------------------------------------------------------


## Weight ------------------------------------------------------------------

long.prey.df %>% 
  filter(Measurement == "Weight_g") %>% 
  filter(Prey_group %in% c("Amphipods","Chironomid","Isopods","Mysid")) %>%
ggplot(aes(x = Species, y = Value, fill = Species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(filter(long.prey.df, Measurement == "Weight_g" & Prey_group %in% c("Amphipods","Chironomid","Isopods","Mysid"))$Value, c(0,0.975))) +
  scale_y_continuous(breaks = seq(0,2,0.25), expand = c(0,0)) +
  scale_fill_jco(labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  facet_wrap(vars(Prey_group), nrow = 1) +
  ylab("Weight (g)") +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )
# ggsave("figures/stomach_contents/major_prey_weight_boxplot_wo_outliers.png", device = "png", dpi = "retina", width = 9.7, height = 4.85, units = "in")


## Count -------------------------------------------------------------------


long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  filter(Prey_group %in% c("Amphipods")) %>% 
ggplot(aes(x = Species, y = Value, fill = Species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(filter(long.prey.df, Measurement == "Count" & Prey_group == "Amphipods")$Value, c(0,0.84))) +
  scale_y_continuous(breaks = seq(0,70,10), expand = c(0,1)) +
  scale_fill_manual(values = cbPalette[c(3,1,2,6)], labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  facet_wrap(vars(Prey_group), nrow = 1) +
  ylab("Count") +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )
# ggsave("figures/stomach_contents/amphipod_count_boxplot_wo_outliers.png", device = "png", dpi = "retina", width = 2.6, height = 4.85, units = "in")

long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  filter(Prey_group %in% c("Chironomid")) %>% 
ggplot(aes(x = Species, y = Value, fill = Species)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(breaks = seq(0,300,50), expand = c(0,1)) +
  scale_fill_manual(values = cbPalette[c(3,1,2,6)], labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  facet_wrap(vars(Prey_group), nrow = 1) +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )
# ggsave("figures/stomach_contents/chironomid_count_boxplot_wo_outliers.png", device = "png", dpi = "retina", width = 2.4, height = 4.85, units = "in")

long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  filter(Prey_group %in% c("Isopods")) %>% 
ggplot(aes(x = Species, y = Value, fill = Species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(filter(long.prey.df, Measurement == "Count" & Prey_group == "Isopods")$Value, c(0,0.92))) +
  scale_y_continuous(breaks = seq(0,40,5), expand = c(0,1)) +
  scale_fill_manual(values = cbPalette[c(3,1,2,6)], labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  facet_wrap(vars(Prey_group), nrow = 1) +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
    )
# ggsave("figures/stomach_contents/isopod_count_boxplot_wo_outliers.png", device = "png", dpi = "retina", width = 2.4, height = 4.85, units = "in")

long.prey.df %>% 
  filter(Measurement == "Count") %>% 
  filter(Prey_group %in% c("Mysid")) %>% 
ggplot(aes(x = Species, y = Value, fill = Species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(filter(long.prey.df, Measurement == "Count" & Prey_group == "Mysid")$Value, c(0,0.85))) +
  scale_y_continuous(breaks = seq(0,40,5), expand = c(0,1)) +
  scale_fill_manual(values = cbPalette[c(3,2,6)], labels = c("Arctic Cisco", "Humpback Whitefish", "Least Cisco")) +
  facet_wrap(vars(Prey_group), nrow = 1) +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
    )
# ggsave("figures/stomach_contents/mysid_count_boxplot_wo_outliers.png", device = "png", dpi = "retina", width = 2.4, height = 4.85, units = "in")


## Relative percent --------------------------------------------------------


long.prey.df %>% 
  filter(Measurement == "Relative_percent") %>% 
  filter(Prey_group %in% c("Amphipods","Mysid","Isopods","Chironomid")) %>% 
ggplot(aes(x = Species, y = Value, fill = Species)) +
  geom_boxplot() +
  scale_fill_manual(values = cbPalette[c(3,1,2,6)], labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10), expand = c(0,1)) +
  facet_wrap(vars(Prey_group), nrow = 1) +
  ylab("Relative Percent") +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )
# ggsave("figures/stomach_contents/major_prey_percent_boxplot_wo_outliers.png", device = "png", dpi = "retina", width = 9.7, height = 4.85, units = "in")

# ANOVA for differences in diet by predator -------------------------------

# glm(predator.species ~ prey.count + prey.weight + prey.percent + ...)

stomach.glm.df <- stomach.df %>% 
  filter(spp_ID != "HBWF_18" & spp_ID != "HBWF_31") %>% 
  filter(Prey_present == 1) %>% 
  dplyr:: select(-c(spp_ID, ID, Prey_present, Relative_fullness, Total_contents_weight_g, Unidentifiable.Count, Vegetation.Count))

stomach.glm.df$Isopods.Count <- as.integer(stomach.glm.df$Isopods.Count)
stomach.glm.df$Amphipods.Count <- as.integer(stomach.glm.df$Amphipods.Count)
stomach.glm.df$Bivalves.Count <- as.integer(stomach.glm.df$Bivalves.Count)
stomach.glm.df$Vertebrate.Weight_g <- as.numeric(stomach.glm.df$Vertebrate.Weight_g)

stomach.glm.df[is.na(stomach.glm.df)] <- 0

levels(stomach.glm.df$Species)
# levels(stomach.glm.df$spp)
# stomach.glm.df$spp <- relevel(stomach.glm.df$Species, ref = "BDWF")


stomach.glm.df$dummy <- ifelse(stomach.glm.df$Species == "ARCS", 1,
                               ifelse(stomach.glm.df$Species == "LSCS", 2,
                                      ifelse(stomach.glm.df$Species == "BDWF", 3, 4)))
test.df <- filter(stomach.glm.df, dummy == 1 | dummy == 2)
test.df$response <- ifelse(test.df$dummy == 1, 1, 0)
attach(test.df)
test.glm <- glm(response ~ . - Species - dummy, family = binomial, data = test.df)
summary(test.glm)
AIC(test.glm)

test2.glm <- glm(response ~ Amphipods.Count, family = binomial, data = test.df)
summary(test2.glm)

test3.glm <- glm(response ~ Isopods.Count, family = binomial, data = test.df)
summary(test3.glm)

AIC(test2.glm, test3.glm, test4.glm, test5.glm)

test4.glm <- glm(response ~ 1, family = binomial, data = test.df)

test5.glm <- glm(response ~ Chironomid.Count, family = binomial, data = test.df)
summary(test5.glm)

arcs.glm.df <- filter(stomach.glm.df, Species == "ARCS")


glm(data = arcs.glm.df)

glm(Species ~ ., data = stomach.glm.df)
#






stomach.glm.df <- droplevels(stomach.glm.df)
test <- multinom(Species ~ 
                   Unidentifiable.Weight_g + Unidentifiable.Relative_percent + Amphipods.Count, data = stomach.glm.df)
summary(test)

z <- summary(test)$coefficients/summary(test)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

exp(coef(test))

head(pp <- fitted(test))
str(stomach.glm.df)
glm.fit <- multinom(Species ~ .-ID, data = stomach.glm.df)
summary(glm.fit)
#Prediction
predict(glm.fit, newdata, "probs")




# NMDS --------------------------------------------------------------------

require(vegan)
require(bugR)

ord1 <- prey.count %>% 
  dplyr::select(-c(Species,ID)) %>% 
  filter(spp_ID != "ARCS_106")
rownames(ord1) <- ord1$spp_ID
ord1 <- ord1 %>% 
  dplyr::select(-spp_ID)
ord1[is.na(ord1)] <- 0
ord1 <- as.matrix(ord1)
ord2 <- ord1[which(rowSums(ord1) > 0),]

weight.ord <- prey.weight %>% 
  dplyr::select(-c(Species,ID,Unidentifiable.Weight_g,Vegetation.Weight_g)) %>% 
  filter(spp_ID != "ARCS_106")
rownames(weight.ord) <- weight.ord$spp_ID
weight.ord <- weight.ord %>% 
  dplyr::select(-spp_ID)
weight.ord[is.na(weight.ord)] <- 0
weight.ord <- as.matrix(weight.ord)
weight.ord <- weight.ord[which(rowSums(weight.ord) > 0),]

env1 <- prey.weight %>% 
  dplyr::select(c(spp_ID, Species)) %>% 
  filter(spp_ID != "ARCS_106")
rownames(env1) <- env1$spp_ID
env1 <- env1 %>% 
  dplyr::select(-spp_ID)

##### Prepare dataset for ordination #####

## Get coefficient of variation in data by rows and columns
cv(ord2)
cv(weight.ord)
## NOTE: Ideally (according to McCune & Grace 2002 textbook), we want <100 for both
## Column cv violates this heavily. Let's deal with this.

## Delete rare species
ord3 <- delRare(ord2)
cv(ord3)

weight.ord2 <- delRare(weight.ord)
cv(weight.ord2)
## NOTE: Default here is getting rid of species not present in >5% of samples
## CV reduced substantially, but still unacceptably high.

## Relativize abundance by column (species)
ord4 <- bugR::rel(ord3)
cv(ord4)

weight.ord3 <- bugR::rel(weight.ord2)
cv(weight.ord3)
## NOTE: This will turn every column into a vector of proportion data, from 0 to 1.
## That gets our CV down some more. Good!



##### Run ordinations #####

## Run a quick and dirty ordination to create a step down plot
NMS(ord4,maxruns=1000,stepdown = TRUE)
## NOTE: This helps identify the ideal dimensionality of the ordination (usually 2D or 3D).
## In these plots, lower stress is better, but we are looking for an inflection point.
## The NMS function creates a folder called "NMS Output" in your working directory where
## it saves the NMS output (matrices of output species and site data).

## Run a full ordination using only 2 (and 3) dimensions
NMS(weight.ord3, maxruns = 1000, stepdown = FALSE, only23 = TRUE)

## Read in site and species point data from ordination output and plot
site1 <- read.csv('NMS Output/NMSPoints2D.csv', row.names = 1)
spp1 <- read.csv('NMS Output/NMSSpecies2D.csv', row.names = 1)
plot(site1)
## NOTE: Great, we have a plot. But what does it mean? Let's dig deeper.



##### Fit environmental data to ordination #####

env1 <- filter(env1, rownames(env1) %in% rownames(site1))
## Compute fit, get vectors
fit1 <- envfit(site1, env1)
## NOTE: Some of our environmental parameters seem not very important to the ordination.

## Get R2 of ordination
axisR2(site1, ord3)


## Rotate ordination axes to align Axis 1 along strongest environmental gradient.
rot1 <- ordRotate(site1, fit1, 'DistancePrimary', x.axis = FALSE, flip = TRUE)
fit2 <- envfit(rot1, env1)
plot(rot1)

## Overlay environmental vectors to create a biplot
for(i in 1:4){
  arrows(0, 0, 0.9 * fit1$factors$centroids[i,1], 0.9 * fit1$factors$centroids[i,2], lty = i + 1)
}
legend('topright', legend = c('ARCS', 'BDWF','HBWF','LSCS'), lty = 2:5,
       bty = 'n')

## Add ellipses for site groupings
col1 <- c('blue', 'green', 'orange', 'red')
lty1 <- c(2:5)
grp1 <- LETTERS[1:4]
ordiellipse(site1, env1$Species, col = col1, lty = lty1)
legend('topright', legend = c('ARCS', 'LSCS', 'BDWF', 'HBWF'), lty = lty1, col = col1, 
       bty = 'n')

## Re-plot for time groupings, with no individual points (for cleanliness)
plot(rot1, type = 'n')
col2 <- c('darkred', 'red', 'orange', 'green', 'darkgreen', 'blue', 'purple')
lty2 <- c(7, 6, 5, 4, 3, 2, 1)
ordiellipse(rot1, env1$Day, col = col2, lty = lty2)
legend('topright', legend = paste0('Day', unique(env1$Day)), lty = lty2, col = col2, 
       bty = 'n')




# NMDS (or MDS?)
# Matrix 1:
# Columns = Group & Gut, where Group is the Predator.Species and Gut is the unique ID for the stomach

# Matrix 2:
# Columns = Gut & all the prey columns (e.g., Amphipod.Count, Isopod.Count, Amphipod.Weight, ...)



# GLM models --------------------------------------------------------------

count.dat <- filter(long.prey.df, Measurement == "Count")
weight.dat <- filter(long.prey.df, Measurement == "Weight_g" & Value > 0)
percent.dat <- filter(long.prey.df, Measurement == "Relative_percent")
percent.dat$Value2 <- percent.dat$Value/100
prey <- c("Amphipods","Chironomid","Isopods","Mysid")

i=1
count.aic.df <- data.frame("n.params" = c(0,1,4,5,9), "amphipods" = rep(NA,5), "chironomid" = rep(NA,5), "isopods" = rep(NA,5), "mysid" = rep(NA,5))
for (i in 1:4) {
  ct.mod1 <- glm.nb(Value ~ 1, data = count.dat, subset = Prey_group == prey[i])
  ct.mod2 <- glm.nb(Value ~ 0 + length_mm, data = count.dat, subset = Prey_group == prey[i])
  ct.mod3 <- glm.nb(Value ~ 0 + Species, data = count.dat, subset = Prey_group == prey[i])
  ct.mod4 <- glm.nb(Value ~ 0 + Species + length_mm, data = count.dat, subset = Prey_group == prey[i])
  ct.mod5 <- glm.nb(Value ~ 0 + Species * length_mm, data = count.dat, subset = Prey_group == prey[i])
  count.aic.df[,i+1] <- AIC(ct.mod1, ct.mod2, ct.mod3, ct.mod4, ct.mod5)[,2]
}
print(count.aic.df)
# No model was deemed as significantly better than the null model
# Model is not working correctly with the "~ 0 + ...". It is still using one of the species as the intercept term.

i=1
weight.aic.df <- data.frame("n.params" = c(0,1,4,5,9), "amphipods" = rep(NA,5), "chironomid" = rep(NA,5), "isopods" = rep(NA,5), "mysid" = rep(NA,5))
for (i in 1:4) {
  wt.mod1 <- glm(Value ~ 1, data = weight.dat, subset = Prey_group == prey[i], family = Gamma(link = log))
  wt.mod2 <- glm(Value ~ 0 + length_mm, data = weight.dat, subset = Prey_group == prey[i], family = Gamma(link = log))
  wt.mod3 <- glm(Value ~ 0 + Species, data = weight.dat, subset = Prey_group == prey[i], family = Gamma(link = log))
  wt.mod4 <- glm(Value ~ 0 + Species + length_mm, data = weight.dat, subset = Prey_group == prey[i], family = Gamma(link = log))
  wt.mod5 <- glm(Value ~ 0 + Species * length_mm, data = weight.dat, subset = Prey_group == prey[i], family = Gamma(link = log))
  weight.aic.df[,i+1] <- AIC(wt.mod1, wt.mod2, wt.mod3, wt.mod4, wt.mod5)[,2]
}
weight.aic.df
# Significant models: Amphipod & Model3, Isopod and Model4

amph3.glm <- glm(Value ~ 0 + Species, data = weight.dat, subset = Prey_group == "Amphipods", family = Gamma(link = log))
summary(amph3.glm)

iso4.glm <- glm(Value ~ 0 + Species + length_mm, data = weight.dat, subset = Prey_group == "Isopods", family = Gamma(link = log))
summary(iso4.glm)

new.df <- data.frame(Species = c("ARCS", "BDWF", "HBWF", "LSCS"), length_mm = 200)
predict.glm(iso4.glm, newdata = new.df, type = "response", se.fit = TRUE)


i=1
percent.aic.df <- data.frame("n.params" = c(0,1,4,5,9), "amphipods" = rep(NA,5), "chironomid" = rep(NA,5), "isopods" = rep(NA,5), "mysid" = rep(NA,5))
for (i in 1:4) {
  pct.mod1 <- glmmTMB(Value2 ~ 1, data = filter(percent.dat, "Prey_group" == prey[i]), family = betabinomial(link = "logit"))
  pct.mod2 <- glmmTMB(Value2 ~ 0 + length_mm, data = filter(percent.dat, "Prey_group" == prey[i]), family = betabinomial(link = "logit"))
  pct.mod3 <- glmmTMB(Value2 ~ 0 + Species, data = filter(percent.dat, "Prey_group" == prey[i]), family = betabinomial(link = "logit"))
  pct.mod4 <- glmmTMB(Value2 ~ 0 + Species + length_mm, data = filter(percent.dat, "Prey_group" == prey[i]), family = betabinomial(link = "logit"))
  pct.mod5 <- glmmTMB(Value2 ~ 0 + Species * length_mm, data = filter(percent.dat, "Prey_group" == prey[i]), family = betabinomial(link = "logit"))
  percent.aic.df[,i+1] <- AIC(pct.mod1, pct.mod2, pct.mod3, pct.mod4, pct.mod5)[,2]
}
percent.aic.df

i=1
percent.aic.df <- data.frame("n.params" = c(0,1,4,5,9), "amphipods" = rep(NA,5), "chironomid" = rep(NA,5), "isopods" = rep(NA,5), "mysid" = rep(NA,5))
for (i in 1:4) {
  pct.mod1 <- glm(Value2 ~ 1, data = percent.dat, subset = Prey_group == prey[i], family = binomial(link = "logit"))
  pct.mod2 <- glm(Value2 ~ 0 + length_mm, data = percent.dat, subset = Prey_group == prey[i], family = binomial(link = "logit"))
  pct.mod3 <- glm(Value2 ~ 0 + Species, data = percent.dat, subset = Prey_group == prey[i], family = binomial(link = "logit"))
  pct.mod4 <- glm(Value2 ~ 0 + Species + length_mm, data = percent.dat, subset = Prey_group == prey[i], family = binomial(link = "logit"))
  pct.mod5 <- glm(Value2 ~ 0 + Species * length_mm, data = percent.dat, subset = Prey_group == prey[i], family = binomial(link = "logit"))
  percent.aic.df[,i+1] <- AIC(pct.mod1, pct.mod2, pct.mod3, pct.mod4, pct.mod5)[,2]
}
percent.aic.df
# Significant models: NONE!


# ggplot of significant models --------------------------------------------

ggplot(data = filter(weight.dat, Prey_group == "Amphipods"), aes(x = Species, y = Value)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0,1.05), breaks = seq(0,1,0.2), expand = c(0,0)) +
  ylab("Amphipod weight (g)") +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )
# ggsave("figures/stomach_contents/amphipod_weight_boxplot_wo_outliers.png", device = "png", dpi = "retina", width = 9.7, height = 4.85, units = "in")

ggplot(data = filter(weight.dat, Prey_group == "Isopods"), aes(x = length_mm, y = Value, shape = Species, color = Species)) +
  geom_point(aes(shape = Species)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(limits = c(-0.1,4.1), breaks = seq(0,4,0.5), expand = c(0,0)) +
  ylab("Isopod weight (g)") +
  xlab("Length (mm)") +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = c(0.1,0.75)
  )
# ggsave("figures/stomach_contents/isopod_weight_by_length.png", device = "png", dpi = "retina", width = 9.7, height = 4.85, units = "in")


# Prey-specific abundance -------------------------------------------------

weight.dat <- filter(long.prey.df, Measurement == "Weight_g" & Value > 0 & 
                       Prey_group %in% c("Mysid","Isopods","Amphipods","Chironomid","Bivalves","Fish","Vertebrate","Insects"))

PSA.df <- weight.dat %>% 
  group_by(Species, Prey_group) %>% 
  summarise(
    "PSA" = sum(Value)
  )

amph <- weight.dat %>% 
  filter(spp_ID %in% filter(weight.dat, Prey_group == "Amphipods")$spp_ID) %>% 
  group_by(Species) %>% 
  summarise("Total" = sum(Value))
amph
mys <- weight.dat %>% 
  filter(spp_ID %in% filter(weight.dat, Prey_group == "Mysid")$spp_ID) %>% 
  group_by(Species) %>% 
  summarise("Total" = sum(Value))
mys
iso <- weight.dat %>% 
  filter(spp_ID %in% filter(weight.dat, Prey_group == "Isopods")$spp_ID) %>% 
  group_by(Species) %>% 
  summarise("Total" = sum(Value))
iso
chiro <- weight.dat %>% 
  filter(spp_ID %in% filter(weight.dat, Prey_group == "Chironomid")$spp_ID) %>% 
  group_by(Species) %>% 
  summarise("Total" = sum(Value))
chiro
bival <- weight.dat %>% 
  filter(spp_ID %in% filter(weight.dat, Prey_group == "Bivalves")$spp_ID) %>% 
  group_by(Species) %>% 
  summarise("Total" = sum(Value))
bival
fish <- weight.dat %>% 
  filter(spp_ID %in% filter(weight.dat, Prey_group == "Fish")$spp_ID) %>% 
  group_by(Species) %>% 
  summarise("Total" = sum(Value))
fish
insect <- weight.dat %>% 
  filter(spp_ID %in% filter(weight.dat, Prey_group == "Insects")$spp_ID) %>% 
  group_by(Species) %>% 
  summarise("Total" = sum(Value))
insect

PSA.df$TotalH <- NA
PSA.df[PSA.df$Prey_group == "Amphipods",]$TotalH <- filter(PSA.df, Prey_group == "Amphipods")$PSA/amph$Total
PSA.df[PSA.df$Prey_group == "Mysid",]$TotalH <- filter(PSA.df, Prey_group == "Mysid")$PSA/mys$Total
PSA.df[PSA.df$Prey_group == "Isopods",]$TotalH <- filter(PSA.df, Prey_group == "Isopods")$PSA/iso$Total
PSA.df[PSA.df$Prey_group == "Chironomid",]$TotalH <- filter(PSA.df, Prey_group == "Chironomid")$PSA/chiro$Total
PSA.df[PSA.df$Prey_group == "Bivalves",]$TotalH <- filter(PSA.df, Prey_group == "Bivalves")$PSA/bival$Total
PSA.df[PSA.df$Prey_group == "Fish",]$TotalH <- filter(PSA.df, Prey_group == "Fish")$PSA/fish$Total
PSA.df[PSA.df$Prey_group == "Insects",]$TotalH <- filter(PSA.df, Prey_group == "Insects")$PSA/insect$Total

PSA.df$PSAh <- PSA.df$TotalH*100


FO.df <- long.prey.df %>% 
  filter(Measurement == "Weight_g" & Prey_group %in% c("Mysid","Isopods","Amphipods","Chironomid","Bivalves","Fish","Vertebrate","Insects")) %>% 
  group_by(Species) %>% 
  summarise(
    "Isopods" = sum(Prey_group == "Isopods"),
    "Amphipods" = sum(Prey_group == "Amphipods"),
    "Bivalves" = sum(Prey_group == "Bivalves"),
    "Insects" = sum(Prey_group == "Insects"),
    "Chironomid" = sum(Prey_group == "Chironomid"),
    "Mysid" = sum(Prey_group == "Mysid"),
    "Fish" = sum(Prey_group == "Fish"),
  )
FO.df

FO.df <- FO.df %>% 
  gather(key = "Prey_group", value = "Value", 2:8, na.rm = T)
FO.df$freq <- ifelse(FO.df$Species == "ARCS", FO.df$Value/46,
                     ifelse(FO.df$Species == "BDWF", FO.df$Value/22,
                            ifelse(FO.df$Species == "HBWF", FO.df$Value/47,
                                   FO.df$Value/56)))
FO.df$freq100 <- FO.df$freq*100

finalPSA.df <- merge(FO.df,PSA.df, by = c("Species", "Prey_group"))

finalPSA.df[finalPSA.df$Prey_group == "Insects",]$Prey_group <- "Emergent 'flies'"
finalPSA.df[finalPSA.df$Prey_group == "Chironomid",]$Prey_group <- "Chironomids"
finalPSA.df[finalPSA.df$Prey_group == "Mysid",]$Prey_group <- "Mysids"

ggplot(finalPSA.df, aes(x = freq100, y = PSAh, label = Prey_group)) +
  geom_point() +
  geom_hline(yintercept = 50, color = "red", lty = 2) +
  geom_vline(xintercept = 50, color = "red", lty = 2) +
  geom_text_repel() +
  facet_wrap(vars(Species)) +
  ylab("Prey-specific abundance") +
  xlab("Frequency of occurrence") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,25), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,25), expand = c(0,0)) +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    panel.spacing = unit(0.3, units = "in"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    plot.margin = unit(c(0,0.2,0.1,0.1), units = "in")
  )
# ggsave("figures/stomach_contents/major_prey_weight_boxplot_wo_outliers.png", device = "png", dpi = "retina", width = 5, height = 4.85, units = "in")



# Histograms
hist(count.dat$Value)
hist(weight.dat$Value)
hist(percent.dat$Value)
