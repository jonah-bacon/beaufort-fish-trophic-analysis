# =====================================================================
# Otolith chronology analysis
# Jonah Bacon
# 28 February 2023
# =====================================================================


# Load libraries ----------------------------------------------------------

library(tidyverse)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

species.names <- c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")
names(species.names) <- c("ARCS","BDWF","HBWF","LSCS")

# Load data ---------------------------------------------------------------

gatt <- read.csv("data/gatt_Beaufort Whitefish Age.Length.csv", header = TRUE)
names(gatt) <- c("species", "ID", "length_mm", "age")
gatt$sampler = "gatt"

bacon <- read.csv("data/data - Otolith Chronology.csv", header = TRUE)
bacon <- bacon %>% select(spp_ID, species, ID, age_estimate_LD)

lw.df <- read.csv("data/data - Sample Data.csv")
lw.df <- lw.df %>% select(spp_ID, species, ID, length_mm)

bacon <- right_join(lw.df, bacon, by = c("spp_ID", "species", "ID"))
bacon <- bacon %>% select(-spp_ID)
names(bacon) <- c("species", "ID", "length_mm", "age")
bacon$sampler = "bacon"

length.age.df <- rbind(gatt, bacon)


# Visualize data ----------------------------------------------------------

ggplot(length.age.df, aes(x = age, y = length_mm)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab("Fork length (mm)") +
  xlab("Age") +
  scale_shape_manual(values=c(21,22,24,25,4), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,28),breaks=seq(0,25,5), expand = c(0,0)) +
  scale_y_continuous(limits=c(0,550),breaks=seq(0,500,100), expand = c(0,0)) +
  facet_wrap(~species, nrow = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )


## Walford plot ------------------------------------------------------------

mean.lg.age.df <- length.age.df %>% 
  group_by(species,age) %>% 
  summarize(mean.lg = mean(length_mm)) %>% 
  filter(!is.na(age))

write.csv(mean.lg.age.df, file = 'data/mean.lg.age.csv') # Manually create mean lg(t+1)
mean.lg.age.df <- read.csv("data/mean.lg.age.csv", header = TRUE)

LVB.df <- mean.lg.age.df %>% 
  group_by(species) %>% 
  summarize(
    K = -log(lm(mean.lg.tplus1 ~ mean.lg)$coef[2]),
    Linf = lm(mean.lg.tplus1 ~ mean.lg)$coef[1]/(1 - lm(mean.lg.tplus1 ~ mean.lg)$coef[2]),
    to = (lm(log(Linf - mean.lg) ~ age)$coef[1] - log(Linf))/K
    )

