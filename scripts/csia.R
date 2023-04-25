# =====================================================================
# Nitrogen Amino Acid Compound-Specific SIA of eye lenses
# Jonah Bacon
# 22 March 2023
# =====================================================================

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggsci)
library(ggrepel)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

species.names <- c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")
names(species.names) <- c("ARCS","BDWF","HBWF","LSCS")

# Load CSIA data ---------------------------------------------------------------

csia.dat <- read.csv("data/data - CSIA Data.csv", header = TRUE, skip = 1)

## Average aa values across three replicate samples:
mean.csia <- csia.dat %>% 
  group_by(sample_ID) %>% 
  summarise(lens_ID = unique(lens_ID), species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            Ala = mean(Ala, na.rm = T), Gly = mean(Gly, na.rm = T), Val = mean(Val, na.rm = T), Leu = mean(Leu, na.rm = T), iLeu = mean(iLeu, na.rm = T),
            nLeu = mean(nLeu, na.rm = T), Thr = mean(Thr, na.rm = T), Ser = mean(Ser, na.rm = T), Asp = mean(Asp, na.rm = T), Pro = mean(Pro, na.rm = T),
            Glu = mean(Glu, na.rm = T), Met = mean(Met, na.rm = T), Phe = mean(Phe, na.rm = T), Caf = mean(Caf, na.rm = T), Lys = mean(Lys, na.rm = T))

# Load lens measurement data ----------------------------------------------

# load("data/lens_measurements.RData")
# 
# lens.df <- lens.df %>% 
#   select(-c(Measurement_pixels, Magnification, Conversion_factor, Measurement_um, photo_ID, Notes, Diameter_um, Radius_um, SIAposition))

# head(lens.df)
# str(lens.df)

# csia.df <- mean.csia %>% 
#   left_join(lens.df, by = c("lens_ID","layer","spp_ID","species","ID"), keep = FALSE)
# csia.df$SIAposition[csia.df$sample_ID == "LS_36-2_07-06"] = 1424.076005/4

# Weighted means estimates of TP ------------------------------------------


## Produce weighted means estimates ((Hayes et al. 1990, Vander Zanden et al. 2013, 
## Glass et al. 2020)
##
## Methods Section 2.5
## https://www.int-res.com/articles/meps2020/641/m641p195.pdf

weighted.means <- csia.dat %>% 
  group_by(CSIA_ID, sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
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

ARCS.weighted.means <- csia.dat %>% 
  filter(species == "ARCS") %>% 
  group_by(sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

BDWF.weighted.means <- csia.dat %>% 
  filter(species == "BDWF") %>% 
  group_by(sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

LSCS.weighted.means <- csia.dat %>% 
  filter(species == "LSCS") %>% 
  group_by(sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

HBWF.weighted.means <- csia.dat %>% 
  filter(species == "HBWF") %>% 
  group_by(sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

combined <- csia.dat %>% 
  group_by(CSIA_ID, sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            TP = ((trophic - source - beta)/TEF) + 1)

## Combine four species TP estimates into one dataframe:
combined.csia.weighted.means.df <- rbind(
  ARCS.weighted.means,
  LSCS.weighted.means,
  BDWF.weighted.means,
  HBWF.weighted.means
)

clean.csia.df <- filter(combined.csia.weighted.means.df, !is.na(TP))
clean.csia.df <- clean.csia.df %>% 
  left_join(csia.dat, by = c("sample_ID","spp_ID","species","ID"), keep = FALSE)

# Phenylalanine across ontogeny -------------------------------------------

ggplot(mean.csia, aes(x = SIAposition, y = Phe, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression("Phenylalanine"~italic(delta)^15*N)) +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  # scale_y_continuous(limits=c(0,5), breaks = 0:5, expand = c(0,0)) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/phenylalanine.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")


# Phenylalanine across ontogeny, by species -------------------------------

## ARCS --------------------------------------------------------------------
ggplot(filter(mean.csia, species == "ARCS"), aes(x = SIAposition, y = Phe, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression("Phenylalanine"~italic(delta)^15*N)) +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  # scale_y_continuous(limits=c(0,5), breaks = 0:5, expand = c(0,0)) +
  # facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/ARCS.phenylalanine.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")


## BDWF --------------------------------------------------------------------
ggplot(filter(mean.csia, species == "BDWF"), aes(x = SIAposition, y = Phe, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(3,2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression("Phenylalanine"~italic(delta)^15*N)) +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  # scale_y_continuous(limits=c(0,5), breaks = 0:5, expand = c(0,0)) +
  # facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/BDWF.phenylalanine.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")


## HBWF --------------------------------------------------------------------
ggplot(filter(mean.csia, species == "HBWF"), aes(x = SIAposition, y = Phe, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression("Phenylalanine"~italic(delta)^15*N)) +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  # scale_y_continuous(limits=c(0,5), breaks = 0:5, expand = c(0,0)) +
  # facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/HBWF.phenylalanine.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")


## LSCS --------------------------------------------------------------------
ggplot(filter(mean.csia, species == "LSCS"), aes(x = SIAposition, y = Phe, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression("Phenylalanine"~italic(delta)^15*N)) +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  # scale_y_continuous(limits=c(0,5), breaks = 0:5, expand = c(0,0)) +
  # facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/LSCS.phenylalanine.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")



# Basic TP vs SIA position ------------------------------------------------

## Weighted means estimates by species, according to methods above:

ARCS.weighted.means <- csia.dat %>% 
  filter(species == "ARCS") %>% 
  group_by(sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            weighted.means.TP = ((trophic - source - beta)/TEF) + 1)

BDWF.weighted.means <- csia.dat %>% 
  filter(species == "BDWF") %>% 
  group_by(sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            weighted.means.TP = ((trophic - source - beta)/TEF) + 1)

LSCS.weighted.means <- csia.dat %>% 
  filter(species == "LSCS") %>% 
  group_by(sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            weighted.means.TP = ((trophic - source - beta)/TEF) + 1)

HBWF.weighted.means <- csia.dat %>% 
  filter(species == "HBWF") %>% 
  group_by(sample_ID) %>% 
  summarise(species = unique(species), ID = unique(ID), spp_ID = unique(spp_ID), SIAposition = unique(SIAposition),
            source.numerator = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.numerator = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denominator = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denominator = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.numerator/source.denominator,
            trophic = trophic.numerator/trophic.denominator,
            weighted.means.TP = ((trophic - source - beta)/TEF) + 1)


## Combine four species TP estimates into one dataframe:
weighted.means.df <- rbind(
  ARCS.weighted.means,
  LSCS.weighted.means,
  BDWF.weighted.means,
  HBWF.weighted.means
)

clean.weighted.means.df <- filter(weighted.means.df, !is.na(weighted.means.TP))

TP.df <- mean.csia %>% 
  mutate(bradley.TP = 1 + ((Glu - Phe - 3.6)/5.7),
         chikaraishi.TP = 1 + ((Glu - Phe - 3.4)/7.6))

TP.df <- clean.weighted.means.df %>% 
  left_join(TP.df, by = c("sample_ID","spp_ID","species","ID","SIAposition"), keep = FALSE) %>% 
  select(c(sample_ID,lens_ID,spp_ID,species,ID,SIAposition,weighted.means.TP, bradley.TP, chikaraishi.TP))



# Weighted means TP vs SIA position ---------------------------------------

# Plot data ---------------------------------------------------------------

# Bradley -----------------------------------------------------------------

ggplot(TP.df, aes(x = SIAposition, y = bradley.TP, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab("Bradley et al. TP") +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(0,5.5), breaks = 0:5, expand = c(0,0)) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/Bradley.CSIA_TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

# Chikaraishi -------------------------------------------------------------

ggplot(TP.df, aes(x = SIAposition, y = chikaraishi.TP, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab("Chikaraishi et al. TP") +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(0,5.5), breaks = 0:5, expand = c(0,0)) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/chikaraishi.CSIA_TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Weighted Means ----------------------------------------------------------

ggplot(TP.df, aes(x = SIAposition, y = weighted.means.TP, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab("Weighted Means TP") +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(0,5.5), breaks = 0:5, expand = c(0,0)) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/weighted.means.CSIA_TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")


# 3 methods per species ---------------------------------------------------

## ARCS
TP.df %>% 
  pivot_longer(cols = c(weighted.means.TP, bradley.TP, chikaraishi.TP), names_to = "method", values_to = "TP_estimate") %>% 
  filter(species == "ARCS") %>% 
ggplot(aes(x = SIAposition, y = TP_estimate, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab("TP estimate") +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(0,5.5), breaks = 0:5, expand = c(0,0)) +
  facet_wrap(~method, nrow = 1, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/ARCS.CSIA_TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## BDWF
TP.df %>% 
  pivot_longer(cols = c(weighted.means.TP, bradley.TP, chikaraishi.TP), names_to = "method", values_to = "TP_estimate") %>% 
  filter(species == "BDWF") %>% 
ggplot(aes(x = SIAposition, y = TP_estimate, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(3,2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab("TP estimate") +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(0,5.5), breaks = 0:5, expand = c(0,0)) +
  facet_wrap(~method, nrow = 1, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/BDWF.CSIA_TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## HBWF
TP.df %>% 
  pivot_longer(cols = c(weighted.means.TP, bradley.TP, chikaraishi.TP), names_to = "method", values_to = "TP_estimate") %>% 
  filter(species == "HBWF") %>% 
ggplot(aes(x = SIAposition, y = TP_estimate, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(2,4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab("TP estimate") +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(0,5.5), breaks = 0:5, expand = c(0,0)) +
  facet_wrap(~method, nrow = 1, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/HBWF.CSIA_TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## LSCS
TP.df %>% 
  pivot_longer(cols = c(weighted.means.TP, bradley.TP, chikaraishi.TP), names_to = "method", values_to = "TP_estimate") %>% 
  filter(species == "LSCS") %>% 
ggplot(aes(x = SIAposition, y = TP_estimate, group = spp_ID, fill = species, shape = species)) + 
  geom_line(linewidth = 0.5, lty = 2) +
  geom_point(aes(color = species, fill = species), color = "black", cex = 3) +
  scale_shape_manual(values=c(25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(4)]) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab("TP estimate") +
  scale_x_continuous(limits=c(0,2100), breaks=seq(0,2000,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(0,5.5), breaks = 0:5, expand = c(0,0)) +
  facet_wrap(~method, nrow = 1, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave(file = "figures/LSCS.CSIA_TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

# d15N + bulk vs lens position, faceted by species -------------------------------

load("data/bulk_sia_data.RData")

ggplot(data = sia.df, aes(x = SIAposition, y = d15N)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 3) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1), sec.axis = sec_axis(~ . /3, name = "Trophic Position")) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
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
# ggsave("figures/d15N.across.lens.by.spp.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Chikaraishi
ggplot(data = sia.df, aes(x = SIAposition, y = d15N)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 3, alpha = 0.2) +
  geom_line(data = TP.df, aes(x = SIAposition, y = chikaraishi.TP*3, group = spp_ID), lwd = 0.5, lty = 2) +
  geom_point(data = TP.df, aes(x = SIAposition, y = chikaraishi.TP*3, group = spp_ID, fill = species, shape = species), color = "black", cex = 3) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1), sec.axis = sec_axis(~ . /3, name = "Chikaraishi TP")) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
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
# ggsave("figures/d15N+chik.TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Bradley
ggplot(data = sia.df, aes(x = SIAposition, y = d15N)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 3, alpha = 0.2) +
  geom_line(data = TP.df, aes(x = SIAposition, y = bradley.TP*3, group = spp_ID), lwd = 0.5, lty = 2) +
  geom_point(data = TP.df, aes(x = SIAposition, y = bradley.TP*3, group = spp_ID, fill = species, shape = species), color = "black", cex = 3) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1), sec.axis = sec_axis(~ . /3, name = "Bradley TP")) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
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
# ggsave("figures/d15N+brad.TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Weighted Means
ggplot(data = sia.df, aes(x = SIAposition, y = d15N)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 3, alpha = 0.2) +
  geom_line(data = TP.df, aes(x = SIAposition, y = weighted.means.TP*3, group = spp_ID), lwd = 0.5, lty = 2) +
  geom_point(data = TP.df, aes(x = SIAposition, y = weighted.means.TP*3, group = spp_ID, fill = species, shape = species), color = "black", cex = 3) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1), sec.axis = sec_axis(~ . /3, name = "Weighted Means TP")) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
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
# ggsave("figures/d15N+wm.TP.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")
