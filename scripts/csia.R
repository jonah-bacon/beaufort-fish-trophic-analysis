# =====================================================================
# Nitrogen amino acid compound-specific SIA from eye lenses
# Jonah Bacon
# 12 October 2022
# =====================================================================


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggsci)
library(ggrepel)

# Load data ---------------------------------------------------------------

csia.dat <- read.csv("data/CSIA_data_1.csv", header = TRUE, skip = 1)

## Remove lab standards from analysis:
temp.dat1 <- filter(csia.dat, CSIA_ID != is.na(CSIA_ID))

## Average aa values across three replicate samples:
temp.dat2 <- temp.dat1 %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            Ala = mean(Ala, na.rm = T), Gly = mean(Gly, na.rm = T), Val = mean(Val, na.rm = T), Leu = mean(Leu, na.rm = T), iLeu = mean(iLeu, na.rm = T),
            nLeu = mean(nLeu, na.rm = T), Thr = mean(Thr, na.rm = T), Ser = mean(Ser, na.rm = T), Asp = mean(Asp, na.rm = T), Pro = mean(Pro, na.rm = T),
            Glu = mean(Glu, na.rm = T), Met = mean(Met, na.rm = T), Phe = mean(Phe, na.rm = T), Caf = mean(Caf, na.rm = T), Lys = mean(Lys, na.rm = T))


# Weighted means estimates of TP ------------------------------------------


## Produce weighted means estimates ((Hayes et al. 1990, Vander Zanden et al. 2013, 
## Glass et al. 2020)
##
## Methods Section 2.5
## https://www.int-res.com/articles/meps2020/641/m641p195.pdf

weighted.means <- temp.dat1 %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            del.15N.Gly = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T)) / sum(1/sd(Gly, na.rm = T)), 
            del.15N.Phe = sum(mean(Phe, na.rm = T)/sd(Phe, na.rm = T)) / sum(1/sd(Phe, na.rm = T)), 
            del.15N.Lys = sum(mean(Lys, na.rm = T)/sd(Lys, na.rm = T)) / sum(1/sd(Lys, na.rm = T)),
            del.15N.Ala = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T)) / sum(1/sd(Ala, na.rm = T)), 
            del.15N.Leu = sum(mean(Leu, na.rm = T)/sd(Leu, na.rm = T)) / sum(1/sd(Leu, na.rm = T)), 
            del.15N.Glu = sum(mean(Glu, na.rm = T)/sd(Glu, na.rm = T)) / sum(1/sd(Glu, na.rm = T)))

## Beta and TEF values from:
## Bradley CJ, Wallsgrove NJ, Choy CA, Drazen JC, Hetherington ED, Hoen DK, Popp BN (2015) 
## Trophic position estimates of marine teleosts using amino acid compound specific isotopic analysis. 
## Limnol Oceanogr Methods 13:476−493 

beta = 3.6
TEF = 5.7


# Weighted means TP estimates by species ----------------------------------
## Weighted means estimates by species, according to methods above:

ARCS.weighted.means <- temp.dat1 %>% 
  filter(Species == "ARCS") %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

BDWF.weighted.means <- temp.dat1 %>% 
  filter(Species == "BDWF") %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

LSCS.weighted.means <- temp.dat1 %>% 
  filter(Species == "LSCS") %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.numerator = sum(mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

HBWF1.weighted.means <- temp.dat1 %>% 
  filter(Species == "HBWF" & CSIA_ID == 8) %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

HBWF2.weighted.means <- temp.dat1 %>% 
  filter(Species == "HBWF" & CSIA_ID == 10) %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.numerator = sum(mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

## Combine four species TP estimates into one dataframe:
combined.csia.weighted.means.df <- rbind(
  ARCS.weighted.means,
  LSCS.weighted.means,
  BDWF.weighted.means,
  HBWF1.weighted.means,
  HBWF2.weighted.means
)

## Create a dummy variable corresponding to the relative time point
## the sample came from within the individual fish's ontogeny
## i.e., 1 = early time point = sample came from the nucleus
## 2 = later time point = sample came from external layer of lens
## Variable values are specific to the order of rbind-ing done when creating
## the "combined.csia.weighted.means.df"
combined.csia.weighted.means.df$tempvar <- c(2,1,2,1,3,1,2,2,1)


## Plot weighted means TP estimates versus relative time point:

ggplot(combined.csia.weighted.means.df, aes(x = tempvar, y = TP, color = Species, fill = Species, label = round(TP,2))) +
  scale_color_jco() +
  scale_fill_jco() +
  geom_line(lwd = 2) +
  geom_point(pch = 21, size = 4, color = "black") +
  geom_label_repel(fill = "white", color = "black", box.padding = 1) +
  xlab("Relative Time Point") +
  ylab("Trophic Position") +
  scale_x_continuous(limits = c(0.8,3.2), breaks = seq(1,3,1), expand = c(0,0)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "gray70"),
    legend.key = element_blank(),
    text = element_text(size = 14, face = "bold")
  )
# ggsave(file = "figures/CSIA_TP.png", width = 7.5, height = 10, units = "in", dpi = 600)


# Load bulk SIA data ------------------------------------------------------


bulk.data <- read.csv("data/bulk.SIA.data.csv", header = T)

delta.df <- bulk.data %>% 
  select(c(species, ID, laminae, d15N, d13C)) %>% 
  unite("spp_ID", species:ID, sep= "_", remove = FALSE) %>% 
  na.omit()

delta.df$spp_ID <- as.factor(delta.df$spp_ID)
delta.df$species <- as.factor(delta.df$species)
delta.df$ID <- as.factor(delta.df$ID)

delta.df %>% 
  filter(spp_ID %in% c("ARCS_14-2", "LSCS_36", "BDWF_38", "HBWF_30-2")) %>% 
  ggplot(aes(x = laminae, y = d15N, color = species, fill = species)) +
  scale_color_jco() +
  scale_fill_jco() +
  geom_line(lwd = 2) +
  geom_point(pch = 21, size = 4, color = "black") +
  # geom_label_repel(fill = "white", color = "black", box.padding = 1) +
  xlab("Lamina layer") +
  ylab("Delta 15 N") +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.2)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "gray70"),
    legend.key = element_blank(),
    text = element_text(size = 14, face = "bold")
  )


# Plot source AA (Phe) vs lamina time -------------------------------------

temp.dat3 <- temp.dat1 %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            Phe_d15N = mean(Phe, na.rm = T),
            Glu_d15N = mean(Glu, na.rm = T))
temp.dat3$tempvar <- c(2,1,2,1,3,1,2,3,1,2)

ggplot(temp.dat3, aes(x = tempvar, y = Phe_d15N, color = Species, fill = Species, label = round(Phe_d15N,2))) +
  scale_color_jco() +
  scale_fill_jco() +
  geom_line(lwd = 2) +
  geom_point(pch = 21, size = 4, color = "black") +
  geom_label_repel(fill = "white", color = "black", box.padding = 1) +
  xlab("Relative Time Point") +
  ylab("Source AA (Phe) d15N") +
  scale_x_continuous(limits = c(0.8,3.2), breaks = seq(1,3,1), expand = c(0,0)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "gray70"),
    legend.key = element_blank(),
    text = element_text(size = 14, face = "bold")
  )

ggplot(temp.dat3, aes(x = tempvar, y = Glu_d15N, color = Species, fill = Species, label = round(Glu_d15N,2))) +
  scale_color_jco() +
  scale_fill_jco() +
  geom_line(lwd = 2) +
  geom_point(pch = 21, size = 4, color = "black") +
  geom_label_repel(fill = "white", color = "black", box.padding = 1) +
  xlab("Relative Time Point") +
  ylab("Trophic AA (Glu) d15N") +
  scale_x_continuous(limits = c(0.8,3.2), breaks = seq(1,3,1), expand = c(0,0)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "gray70"),
    legend.key = element_blank(),
    text = element_text(size = 14, face = "bold")
  )
