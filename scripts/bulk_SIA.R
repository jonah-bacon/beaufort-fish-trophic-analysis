# =====================================================================
# Bulk Carbon & Nitrogen SIA of eye lenses
# Jonah Bacon
# 20 May 2022
# =====================================================================

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggsci)
library(nlme)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

species.names <- c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")
names(species.names) <- c("ARCS","BDWF","HBWF","LSCS")

# Reference standards -----------------------------------------------------

stds.bulk.df <- read.csv("data/bulk.SIA.standards.csv", header = T)
stds.bulk.df <- select(stds.bulk.df, -"Notes")
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


# Load bulk SIA data ------------------------------------------------------

bulk.df <- read.csv("data/bulk.SIA.data.csv", header = T)

bulk.df <- bulk.df %>% 
  select(-c(laminae, data_sheet_ID, well_ID, tray_ID, Sample.Name, Sample.Wt., N.Signal, C.Signal, Conc.C, Conc.N, Notes)) %>% 
  unite("image_ID", species:layer, sep = "_", remove = FALSE) %>% 
  unite("lens_ID", species:ID, sep= "_", remove = FALSE)

# head(bulk.df)
# str(bulk.df)

# Load lens measurement data ----------------------------------------------

load("data/lens_measurements.RData")

lens.df <- lens.df %>% 
  select(-c(Measurement_pixels, Magnification, Conversion_factor, Measurement_um))

# head(lens.df)
# str(lens.df)

# Merge data frames -------------------------------------------------------

sia.df <- bulk.df %>% 
  full_join(lens.df, by = c("lens_ID", "layer"), keep = FALSE) %>% 
  filter(!is.na(d15N) | !is.na(d13C))

sia.df$lens_ID <- as.factor(sia.df$lens_ID)
sia.df$spp_ID <- as.factor(sia.df$spp_ID)
sia.df$species <- as.factor(sia.df$species)
sia.df$ID <- as.factor(sia.df$ID)

arcs.bulk.df <- sia.df %>% filter(species == "ARCS")
lscs.bulk.df <- sia.df %>% filter(species == "LSCS")
bdwf.bulk.df <- sia.df %>% filter(species == "BDWF")
hbwf.bulk.df <- sia.df %>% filter(species == "HBWF")

# Visualize the data ------------------------------------------------------

## d15N vs lens position, faceted by species -------------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = d15N)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 3) +
  xlab("Lens radial position (um)") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0.1,3600), breaks=seq(0,4000,500), expand=c(0,5)) +
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0.1)) +
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
# ggsave("figures/d15N.across.lens.by.spp.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## d13C vs lens position, faceted by species -------------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = d13C)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 3) +
  xlab("Lens radial position (um)") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0.1,3600), breaks=seq(0,4000,500), expand=c(0,5)) +
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
# ggsave("figures/d15N.across.lens.by.spp.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Trend in d15N across lens position --------------------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = d15N, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  xlab("Lens radial position (um)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0.1,3600),breaks=seq(0,3600,1000), expand = c(0,5)) +
  scale_y_continuous(limits=c(0,15),breaks=seq(0,15,3), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.85,0.1),
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    axis.line=element_line()
  )
# ggsave("figures/d15N.across.lens.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(data = sia.df, aes(x = SIAposition, y = d15N, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5, alpha = 0.2, show.legend = F) +
  geom_smooth(size = 2, method = "nls", formula = y ~ a - b*exp(-c*x), method.args = list(start = c(a = 13.7, b = 7.514, c= 0.0006322)), se = FALSE) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  xlab("Lens radial position (um)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_x_continuous(limits=c(0.1,3600),breaks=seq(0,3600,1000), expand = c(0,5)) +
  scale_y_continuous(limits=c(0,15),breaks=seq(0,15,3), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.85,0.1),
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    axis.line=element_line()
  )
# ggsave("figures/d15N.across.lens.wtrend.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Trend in d13C across lens position --------------------------------------

ggplot(data = sia.df, aes(x = SIAposition, y = d13C, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  xlab("Lens radial position (um)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(0.1,3600),breaks=seq(0,3600,1000), expand = c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.85,0.1),
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    axis.line=element_line()
  )
# ggsave("figures/d13C.across.lens.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(data = sia.df, aes(x = SIAposition, y = d13C, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5, alpha = 0.2, show.legend = F) +
  geom_smooth(size = 2, method = "nls", formula = y ~ a - b*exp(-c*x), method.args = list(start = c(a = -20, b = 12, c= 0.0005)), se = FALSE) +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  xlab("Lens radial position (um)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_x_continuous(limits=c(0.1,3600),breaks=seq(0,3600,1000), expand = c(0,5)) +
  scale_y_continuous(limits=c(-35,-19),breaks=seq(-35,-20,5), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.85,0.1),
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    axis.line=element_line()
  )
# ggsave("figures/d13C.across.lens.wtrend.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")



# Likelihood ratio test of models -----------------------------------------

# Filter data to remove instances where SIAposition = 0
nls.df <- filter(sia.df, SIAposition != 0)

## For d15N
# Fit nonlinear functions to data, one for each species
model.arcs <- nls(d15N ~ a1 - b1 * exp(-c1 * SIAposition), data = subset(nls.df, species == "ARCS"), start = list(a1 = 15, b1 = 3, c1 = 0.0005))
model.bdwf <- nls(d15N ~ a2 - b2 * exp(-c2 * SIAposition), data = subset(nls.df, species == "BDWF"), start = list(a2 = 12, b2 = 3, c2 = 0.0005))
model.hbwf <- nls(d15N ~ a3 - b3 * exp(-c3 * SIAposition), data = subset(nls.df, species == "HBWF"), start = list(a3 = 15, b3 = 3, c3 = 0.0005))
model.lscs <- nls(d15N ~ a4 - b4 * exp(-c4 * SIAposition), data = subset(nls.df, species == "LSCS"), start = list(a4 = 15, b4 = 3, c4 = 0.0005))

# Perform likelihood ratio test
lrtest <- anova(model.arcs, model.bdwf, model.hbwf, model.lscs)

# Print likelihood ratio test results
print(lrtest)

## For d13C
# Fit nonlinear functions to data, one for each species
model.arcs <- nls(d13C ~ a1 - b1 * exp(-c1 * SIAposition), data = subset(nls.df, species == "ARCS"), start = list(a1 = -22, b1 = 5, c1 = 0.0005))
model.bdwf <- nls(d13C ~ a2 - b2 * exp(-c2 * SIAposition), data = subset(nls.df, species == "BDWF"), start = list(a2 = -22, b2 = 5, c2 = 0.0005))
model.hbwf <- nls(d13C ~ a3 - b3 * exp(-c3 * SIAposition), data = subset(nls.df, species == "HBWF"), start = list(a3 = -22, b3 = 5, c3 = 0.0005))
model.lscs <- nls(d13C ~ a4 - b4 * exp(-c4 * SIAposition), data = subset(nls.df, species == "LSCS"), start = list(a4 = -22, b4 = 5, c4 = 0.0005))

# Perform likelihood ratio test
lrtest <- anova(model.arcs, model.bdwf, model.hbwf, model.lscs)

# Print likelihood ratio test results
print(lrtest)

