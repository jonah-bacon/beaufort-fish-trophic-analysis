# =====================================================================
# Bulk Carbon & Nitrogen SIA of eye lenses
# Jonah Bacon
# 27 February 2023
# =====================================================================

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggsci)
library(nlme)
library(broom)
library(stats)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

species.names <- c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")
names(species.names) <- c("ARCS","BDWF","HBWF","LSCS")

# Load bulk SIA data ------------------------------------------------------

bulk.df <- read.csv("data/data - Bulk SIA Data.csv", header = T)

bulk.df <- bulk.df %>% 
  filter(!Sample.Name %in% c("F28","C30","A32","C31","B37","B76","B77","B42","B43")) %>% # Remove 9 samples where the photo malfunctioned
  filter(!Sample.Name %in% c("C58")) %>% # Remove 1 sample too small to obtain reliable data
  dplyr::select(-c(laminae, data_sheet_ID, well_ID, tray_ID, Sample.Name, Sample.Wt., N.Signal, C.Signal, Notes))

# head(bulk.df)
# str(bulk.df)

# Load lens measurement data ----------------------------------------------

load("data/lens_measurements.RData")

lens.df <- lens.df %>% 
  dplyr::select(-c(Measurement_pixels, Magnification, Conversion_factor, Measurement_um, photo_ID, Notes, Diameter_um, Radius_um))

# head(lens.df)
# str(lens.df)

# Merge data frames -------------------------------------------------------

sia.df <- bulk.df %>% 
  left_join(lens.df, by = c("lens_ID","layer","spp_ID","species","ID"), keep = FALSE) %>% 
  filter(!is.na(d15N) | !is.na(d13C))

sia.df$lens_ID <- as.factor(sia.df$lens_ID)
sia.df$spp_ID <- as.factor(sia.df$spp_ID)
sia.df$species <- as.factor(sia.df$species)
sia.df$ID <- as.factor(sia.df$ID)

# arcs.bulk.df <- sia.df %>% filter(species == "ARCS")
# lscs.bulk.df <- sia.df %>% filter(species == "LSCS")
# bdwf.bulk.df <- sia.df %>% filter(species == "BDWF")
# hbwf.bulk.df <- sia.df %>% filter(species == "HBWF")

save.image(file = "data/bulk_sia_data.RData")
# load("data/bulk_sia_data.RData")

# Reference standards -----------------------------------------------------

stds.bulk.df <- read.csv("data/bulk.SIA.standards.csv", header = T)
stds.bulk.df <- dplyr::select(stds.bulk.df, -Notes)
stds.bulk.df <- na.omit(stds.bulk.df)

stds.summary_total <- data.frame(
  metric = c("mean", "std dev", "expected", "difference"),
  concN = c(mean(stds.bulk.df$Conc.N), sd(stds.bulk.df$Conc.N), 15.30, mean(stds.bulk.df$Conc.N) - 15.30),
  concC = c(mean(stds.bulk.df$Conc.C), sd(stds.bulk.df$Conc.C), 44.30, mean(stds.bulk.df$Conc.C) - 44.30),
  d15N = c(mean(stds.bulk.df$d15N), sd(stds.bulk.df$d15N), 7.00, mean(stds.bulk.df$d15N) - 7.00),
  d13C = c(mean(stds.bulk.df$d13C), sd(stds.bulk.df$d13C), -15.80, mean(stds.bulk.df$d13C) + 15.80))

stds.summary_total[1,] # Observed average value for reference check standards
stds.summary_total[2,] # Standard deviation between reference check standards
stds.summary_total[3,] # Expected value for reference check standards
stds.summary_total[4,] # Difference between observed and expected values for reference check standards

# 95% CI of d15N and d13C values
stds.summary_total[2,2:5]*1.96
# Summary: d15N values of my samples are 'truly' different from each other if they are >0.5261656 (2*1.96*sd) per mille different
# Summary: d13C values of my samples are 'truly' different from each other if they are >0.2799874 (2*1.96*sd) per mille different


# Analysis of standards by each tray of samples that were run:

# stds.summary_mean.by.run <- stds.bulk.df %>% 
#   group_by(tray_ID) %>% 
#   summarise(
#     metric = "mean",
#     concN = mean(Conc.N), 
#     concC = mean(Conc.C), 
#     d15N = mean(d15N), 
#     d13C = mean(d13C))

# stds.summary_sdev.by.run <- stds.bulk.df %>% 
#   group_by(tray_ID) %>% 
#   summarise(
#     metric = "std dev",
#     concN = sd(Conc.N),
#     concC = sd(Conc.C),
#     d15N = sd(d15N),
#     d13C = sd(d13C))

# stds.summary_diff.by.run <- stds.bulk.df %>% 
#   group_by(tray_ID) %>% 
#   summarise(
#     metric = "difference from expected",
#     concN = mean(Conc.N) - 15.30,
#     concC = mean(Conc.C) - 44.30,
#     d15N = mean(d15N) - 7.00,
#     d13C = mean(d13C) + 15.80)

# Visualize the data ------------------------------------------------------

## Distribution of tissue sampling variance -----------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = SIArange)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 2) +
  xlab("Sample lens axis position") +
  ylab("Sample position range") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,1000), breaks = seq(0,1000, 200), expand = c(0,5)) +
  facet_wrap(~species, ncol = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave("figures/lens.sampling.variance.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

var(sia.df$SIArange, na.rm = T)/var(sia.df$SIAposition, na.rm = T)

## C:N ratios --------------------------------------------------------------

ggplot(data = sia.df, aes(Conc.C/Conc.N)) +
  geom_histogram(binwidth = 0.05, color = "white") + 
  geom_vline(aes(xintercept = 3.5), color = "red", linewidth = 1, lty = 2) +
  ylab("Number of samples") +
  xlab("C:N ratio") +
  scale_y_continuous(limits = c(0,110), breaks = seq(0,100,20), expand = c(0,0)) +
  scale_x_continuous(limits = c(2.3,3.6), breaks = seq(2.50, 3.50, 0.25), expand = c(0,0)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    axis.line=element_line()
  )
# ggsave("figures/CN.ratio.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## d15N vs lens position, faceted by species -------------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = d15N)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 3) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1)) +
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

## d13C vs lens position, faceted by species -------------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = d13C)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 3) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(-35,-20), breaks=seq(-35,-20, 5), expand=c(0,0.1)) +
  facet_wrap(~species, ncol = 1, labeller = labeller(species = species.names)) +
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
# ggsave("figures/d13C.across.lens.by.spp.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Trend in d15N across lens position --------------------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = d15N, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800),breaks=seq(0,1800,500), expand = c(0,5)) +
  scale_y_continuous(limits=c(0,15),breaks=seq(0,15,3), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.85,0.2),
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    axis.line=element_line()
  )
# ggsave("figures/d15N.across.lens.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(data = sia.df, aes(x = SIAposition, y = d15N, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5, alpha = 0.2, show.legend = F) +
  geom_smooth(linewidth = 2, method = "nls", formula = y ~ a - b*exp(-c*x), method.args = list(start = c(a = 13.7, b = 7.514, c= 0.0006322)), se = FALSE) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_x_continuous(limits=c(0,1800),breaks=seq(0,1800,500), expand = c(0,5)) +
  scale_y_continuous(limits=c(0,15),breaks=seq(0,15,3), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.85,0.2),
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    axis.line=element_line()
  )
# ggsave("figures/d15N.across.lens.wtrend.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Trend in d13C across lens position --------------------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = d13C, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800),breaks=seq(0,1800,500), expand = c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.85,0.2),
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    axis.line=element_line()
  )
# ggsave("figures/d13C.across.lens.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(data = sia.df, aes(x = SIAposition, y = d13C, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5, alpha = 0.2, show.legend = F) +
  geom_smooth(size = 2, method = "nls", formula = y ~ a - b*exp(-c*x), method.args = list(start = c(a = -20, b = 12, c= 0.0005)), se = FALSE) +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_x_continuous(limits=c(0,1800),breaks=seq(0,1800,500), expand = c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.85,0.2),
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    axis.line=element_line()
  )
# ggsave("figures/d13C.across.lens.wtrend.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")



## d15N individual tracks --------------------------------------------------

sia.df %>% 
  group_by(lens_ID) %>% 
  filter(length(lens_ID) > 2) %>% 
ggplot(aes(x = SIAposition, y = d15N, color = species, group = spp_ID)) +
  geom_line(aes(color = species), alpha = 0.3) +
  geom_point(aes(fill = NULL, shape = species, color = species), cex = 1, alpha = 0.5) +
  # geom_smooth(linewidth = 2, method = "nls", formula = y ~ a - b*exp(-c*x), method.args = list(start = c(a = 15, b = 15, c= 0.01)), se = FALSE) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1)) +
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
# ggsave("figures/d15N.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## d13C individual tracks --------------------------------------------------

sia.df %>% 
  group_by(lens_ID) %>% 
  filter(length(lens_ID) > 2) %>% 
ggplot(aes(x = SIAposition, y = d13C, color = species, group = spp_ID)) +
  geom_line(aes(color = species), alpha = 0.3) +
  geom_point(aes(fill = NULL, shape = species, color = species), cex = 1, alpha = 0.5) +
  # geom_smooth(linewidth = 2, method = "nls", formula = y ~ a - b*exp(-c*x), method.args = list(start = c(a = 15, b = 15, c= 0.01)), se = FALSE) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
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
# ggsave("figures/d13C.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")


## d15N individual tracks, by species --------------------------------------

## Arctic Cisco
sia.df %>% 
  group_by(lens_ID) %>% 
  filter(species == "ARCS" & length(unique(layer)) > 2) %>% 
ggplot(aes(x = SIAposition, y = d15N)) +
  geom_line(aes(color = species), color = "gray85") +
  geom_point(color = "black", cex = 1) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  # scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  # scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1)) +
  facet_wrap(~spp_ID, nrow = 2) +
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
# ggsave("figures/ARCS.d15N.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Broad Whitefish
sia.df %>% 
  group_by(lens_ID) %>% 
  filter(species == "BDWF" & length(unique(layer)) > 2) %>% 
ggplot(aes(x = SIAposition, y = d15N)) +
  geom_line(aes(color = species), color = "gray85") +
  geom_point(color = "black", cex = 1) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  # scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  # scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1)) +
  facet_wrap(~spp_ID, nrow = 2) +
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
# ggsave("figures/BDWF.d15N.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Humpback Whitefish
sia.df %>% 
  group_by(lens_ID) %>% 
  filter(species == "HBWF" & length(unique(layer)) > 2) %>% 
ggplot(aes(x = SIAposition, y = d15N)) +
  geom_line(aes(color = species), color = "gray85") +
  geom_point(color = "black", cex = 1) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  # scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  # scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1)) +
  facet_wrap(~spp_ID, nrow = 2) +
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
# ggsave("figures/HBWF.d15N.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Least Cisco
sia.df %>% 
  group_by(lens_ID) %>% 
  filter(species == "LSCS" & length(unique(layer)) > 2) %>% 
ggplot(aes(x = SIAposition, y = d15N)) +
  geom_line(aes(color = species), color = "gray85") +
  geom_point(color = "black", cex = 1) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  # scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  # scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1)) +
  facet_wrap(~spp_ID, nrow = 2) +
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
# ggsave("figures/LSCS.d15N.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## d13C individual tracks, by species --------------------------------------

## Arctic Cisco
sia.df %>% 
  group_by(lens_ID) %>% 
  filter(species == "ARCS" & length(unique(layer)) > 2) %>% 
ggplot(aes(x = SIAposition, y = d13C)) +
  geom_line(aes(color = species), color = "gray85") +
  geom_point(color = "black", cex = 1) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  # scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  # scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
  facet_wrap(~spp_ID, nrow = 2) +
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
# ggsave("figures/ARCS.d13C.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Broad Whitefish
sia.df %>% 
  group_by(lens_ID) %>% 
  filter(species == "BDWF" & length(unique(layer)) > 2) %>% 
ggplot(aes(x = SIAposition, y = d13C)) +
  geom_line(aes(color = species), color = "gray85") +
  geom_point(color = "black", cex = 1) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  # scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  # scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
  facet_wrap(~spp_ID, nrow = 2) +
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
# ggsave("figures/BDWF.d13C.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Humpback Whitefish
sia.df %>% 
  group_by(lens_ID) %>% 
  filter(species == "HBWF" & length(unique(layer)) > 2) %>% 
ggplot(aes(x = SIAposition, y = d13C)) +
  geom_line(aes(color = species), color = "gray85") +
  geom_point(color = "black", cex = 1) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  # scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  # scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
  facet_wrap(~spp_ID, nrow = 2) +
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
# ggsave("figures/HBWF.d13C.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Least Cisco
sia.df %>% 
  group_by(lens_ID) %>% 
  filter(species == "LSCS" & length(unique(layer)) > 2) %>% 
ggplot(aes(x = SIAposition, y = d13C)) +
  geom_line(aes(color = species), color = "gray85") +
  geom_point(color = "black", cex = 1) +
  xlab("Radial midpoint of lamina (\u03bcm)") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  # scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  # scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0,1800), breaks=seq(0,1800,500), expand=c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
  facet_wrap(~spp_ID, nrow = 2) +
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
# ggsave("figures/LSCS.d13C.individual.tracks.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")



# Models ------------------------------------------------------------------


# Filter data to remove instances where SIAposition = 0
nls.df <- filter(sia.df, SIAposition != 0)
fish <- c("ARCS","HBWF","LSCS")

i=1
d15N.aic.df <- data.frame("n.params" = c(0,2,3), "ARCS" = rep(NA,3), "HBWF" = rep(NA,3), "LSCS" = rep(NA,3))
for (i in 1:3) {
  d15N.mod1 <- glm(d15N ~ 1, data = nls.df, subset = species == fish[i])
  d15N.mod2 <- glm(d15N ~ SIAposition, data = nls.df, subset = species == fish[i])
  d15N.mod3 <- nls(d15N ~ a1 - b1 * exp(-c1 * SIAposition), data = nls.df, subset = species == fish[i], start = list(a1 = 15, b1 = 3, c1 = 0.0005))
  d15N.aic.df[,i+1] <- AIC(d15N.mod1, d15N.mod2, d15N.mod3)[,2]
}

d15N.mod1 <- glm(d15N ~ 1, data = nls.df, subset = species == "BDWF")
d15N.mod2 <- glm(d15N ~ SIAposition, data = nls.df, subset = species == "BDWF")
d15N.mod3 <- nls(d15N ~ a1 - b1 * exp(-c1 * SIAposition), data = nls.df, subset = species == "BDWF", start = list(a1 = 20, b1 = 10, c1 = 0.0005))
d15N.aic.df$BDWF <- AIC(d15N.mod1, d15N.mod2, d15N.mod3)[,2]
print(d15N.aic.df)

## For d15N
# Fit nonlinear functions to data, one for each species
d15N.model.arcs <- glm(d15N ~ SIAposition, data = subset(nls.df, species == "ARCS"))
d15N.model.bdwf <- glm(d15N ~ SIAposition, data = subset(nls.df, species == "BDWF"))
d15N.model.hbwf <- nls(d15N ~ a3 - b3 * exp(-c3 * SIAposition), data = subset(nls.df, species == "HBWF"), start = list(a3 = 15, b3 = 3, c3 = 0.0005))
d15N.model.lscs <- nls(d15N ~ a4 - b4 * exp(-c4 * SIAposition), data = subset(nls.df, species == "LSCS"), start = list(a4 = 15, b4 = 3, c4 = 0.0005))

# Print estimates and 95% CI of coefficients
coef(d15N.model.arcs)
confint(d15N.model.arcs)

coef(d15N.model.bdwf)
confint(d15N.model.bdwf)

coef(d15N.model.hbwf)
confint(d15N.model.hbwf)

coef(d15N.model.lscs)
confint(d15N.model.lscs)


## For d13C
fish <- c("ARCS","BDWF","HBWF","LSCS")

i=1
d13C.aic.df <- data.frame("n.params" = c(0,2,3), "ARCS" = rep(NA,3), "BDWF" = rep(NA,3), "HBWF" = rep(NA,3), "LSCS" = rep(NA,3))
for (i in 1:4) {
  d13C.mod1 <- glm(d13C ~ 1, data = nls.df, subset = species == fish[i])
  d13C.mod2 <- glm(d13C ~ SIAposition, data = nls.df, subset = species == fish[i])
  d13C.mod3 <- nls(d13C ~ a1 - b1 * exp(-c1 * SIAposition), data = nls.df, subset = species == fish[i], start = list(a1 = 15, b1 = 3, c1 = 0.001))
  d13C.aic.df[,i+1] <- AIC(d13C.mod1, d13C.mod2, d13C.mod3)[,2]
}
print(d13C.aic.df)

# Fit functions to data, one for each species
d13C.model.arcs <- nls(d13C ~ a1 - b1 * exp(-c1 * SIAposition), data = subset(nls.df, species == "ARCS"), start = list(a1 = -22, b1 = 5, c1 = 0.001))
d13C.model.bdwf <- glm(d13C ~ SIAposition, data = subset(nls.df, species == "BDWF"))
d13C.model.hbwf <- glm(d13C ~ SIAposition, data = subset(nls.df, species == "HBWF"))
d13C.model.lscs <- nls(d13C ~ a4 - b4 * exp(-c4 * SIAposition), data = subset(nls.df, species == "LSCS"), start = list(a4 = -22, b4 = 5, c4 = 0.001))

# Print estimates and 95% CI of coefficients
coef(d13C.model.arcs)
summary(d13C.model.arcs)
coef(d13C.model.arcs)[1] + 1.96*1.621712
coef(d13C.model.arcs)[1] - 1.96*1.621712
coef(d13C.model.arcs)[2] + 1.96*1.137070
coef(d13C.model.arcs)[2] - 1.96*1.137070
coef(d13C.model.arcs)[3] + 1.96*0.000513
coef(d13C.model.arcs)[3] - 1.96*0.000513

coef(d13C.model.bdwf)
confint(d13C.model.bdwf)

coef(d13C.model.hbwf)
confint(d13C.model.hbwf)

coef(d13C.model.lscs)
confint(d13C.model.lscs)
