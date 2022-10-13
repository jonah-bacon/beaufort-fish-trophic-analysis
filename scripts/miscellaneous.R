library(tidyverse)

temp.df <- read.csv("data/sample.LW.data.csv", header = T)
head(temp.df)
str(temp.df)
temp.df$species <- as.factor(temp.df$species)

ggplot(temp.df, aes(x = length_mm, fill = species, color = species)) +
  geom_histogram() +
  facet_wrap(vars(species))

temp.df %>% 
  group_by(species) %>% 
  summarise("N.less.120" = length(length_mm[length_mm < 120]), 
            "N.120-249" = length(length_mm[length_mm >= 120 & length_mm < 250]),
            "N.greater.250" = length(length_mm[length_mm >= 250]),
            "total" = length(length_mm))
