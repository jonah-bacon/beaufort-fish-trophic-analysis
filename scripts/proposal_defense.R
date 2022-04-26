# =====================================================================
# Bulk carbon & nitrogen SIA from eye lenses
# Jonah Bacon
# 25 April 2022
# =====================================================================

# Load packages -----------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(ggsci)
# library(ggrepel)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Load eye lens data ------------------------------------------------------

bulk.data <- read.csv("data/bulk.SIA.data.csv", header = T)

head(bulk.data)
tail(bulk.data)
str(bulk.data)

delta.df <- bulk.data %>% 
  select(c(species, ID, laminae, d15N, d13C)) %>% 
  unite("spp_ID", species:ID, sep= "_", remove = FALSE)

delta.df$spp_ID <- as.factor(delta.df$spp_ID)
delta.df$species <- as.factor(delta.df$species)
delta.df$ID <- as.factor(delta.df$ID)

head(delta.df)
tail(delta.df)
str(delta.df)

arcs.delta.df <- delta.df %>% filter(species == "ARCS")
lscs.delta.df <- delta.df %>% filter(species == "LSCS")
bdwf.delta.df <- delta.df %>% filter(species == "BDWF")
hbwf.delta.df <- delta.df %>% filter(species == "HBWF")

species.names <- c("Arctic cisco", "Least cisco", "Broad whitefish", "Humpback whitefish")
names(species.names) <- c("ARCS", "LSCS", "BDWF", "HBWF")


# Load muscle/fin data ----------------------------------------------------

muscle.fin.df <- read.csv("data/wooller Isotope class release 22Mar31.csv", header = F, skip = 9, nrows = 13)

muscle.fin.df <- muscle.fin.df %>% 
  select(-c("V9"))
colnames(muscle.fin.df) <- c("Sample.Name","Sample.wt.mg","N.sig.V","C.sig.V","Conc.N","Conc.C","d15N","d13C")

# muscle.fin.df <- muscle.fin.df %>% 
#   mutate("C:N" = Conc.C/Conc.N)
# 
# standards.df <- muscle.fin.df %>% 
#   filter(Sample.Name == "ref/chk/peptone")

muscle.fin.df <- muscle.fin.df %>% 
  filter(Sample.Name != "ref/chk/peptone") %>% 
  separate("Sample.Name", into = c("Sample.Name", "Sample.Type"), sep = -1, remove = T)

fin.df <- muscle.fin.df %>% 
  filter(Sample.Type == "F") %>% 
  summarise("Sample.Name" = Sample.Name, "fin.d15N" = d15N, "fin.d13C" = d13C)

muscle.df <- muscle.fin.df %>% 
  filter(Sample.Type == "M") %>% 
  summarise("Sample.Name" = Sample.Name, "muscle.d15N" = d15N, "muscle.d13C" = d13C)

muscle.fin.df <- merge(muscle.df, fin.df)

correlation.df <- merge(muscle.fin.df, 
                     delta.df %>% 
                       select(spp_ID, laminae, d15N, d13C) %>% 
                       filter(spp_ID %in% c("BDWF_41", "LSCS_42", "HBWF_33", "LSCS_44", "BDWF_45_1")) %>% 
                       mutate(Sample.Name = spp_ID, Sample.Type = laminae, d15N = d15N, d13C = d13C, .keep = "none") %>% 
                       group_by(Sample.Name) %>% 
                       summarise("lamina.d15N" = d15N[1], "lamina.d13C" = d13C[1], "mean.lamina.d15N" = mean(d15N[1:3]), "mean.lamina.d13C" = mean(d13C[1:3]))
                     )

# Visualize the data ------------------------------------------------------

## Isotope value reconstructions -------------------------------------------

### d13C vs lamina layer, faceted by species --------------------------------

ggplot(data = delta.df, aes(x = laminae, y = d13C)) +
  geom_line(aes(col = ID), cex = 1.2, alpha = 0.5) +
  geom_point(aes(col = ID), cex = 4) +
  xlab("Lamina layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_color_hue() +
  scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1), expand=c(0,0.1)) +
  scale_y_continuous(limits=c(-34,-20), breaks=seq(-34,-20, 2), expand=c(0,0.1)) +
  facet_wrap(~species, ncol = 1, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size=18),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    strip.text = element_text(size = 13),
    axis.line=element_line(),
    legend.position = "none"
  )


### d15N vs lamina layer, faceted by species --------------------------------

ggplot(data = delta.df, aes(x = laminae, y = d15N)) +
  geom_line(aes(col = ID), cex = 1.2, alpha = 0.5) +
  geom_point(aes(col = ID), cex = 4) +
  xlab("Lamina layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_color_hue() +
  scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1), expand=c(0,0.1)) +
  scale_y_continuous(limits = c(3.5,13), breaks = seq(4,13,1), expand = c(0,0.1)) +
  facet_wrap(~species, ncol = 1, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size=18),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    strip.text = element_text(size = 13),
    axis.line=element_line(),
    legend.position = "none"
  )



d15N.df <- fish.df %>%
  select(Sample.Name, Sample.Type, d15N) %>% 
  pivot_wider(names_from = Sample.Type, values_from = d15N)

d15N.intercept.hat <- mean(d15N.df$F - d15N.df$M)
d15N.intercept.hat

d15N.lm <- lm(F ~ M, data = d15N.df)
summary(d15N.lm)

ggplot(d15N.df, aes(x = Muscle, y = Fin, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  geom_abline(slope = 1, intercept = d15N.intercept.hat, cex = 1.5, color = "black") +
  geom_point(shape = 17, size = 5) +
  xlab(expression(Muscle~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  ylab(expression(Fin~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  scale_x_continuous(limits = c(8,14), breaks = seq(8,14,0.5), expand = c(0,0.1)) +
  scale_y_continuous(limits = c(8,14), breaks = seq(8,14,0.5), expand = c(0,0.1)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2")) +
  labs(color = "Sample Name") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    legend.position = c(0.1,0.7),
    legend.background = element_rect(fill = "gray90", color = "black"),
    legend.key = element_rect(fill = "gray90"),
    legend.title = element_text(size =14),
    legend.text = element_text(size=13),
    axis.line=element_line()
  )


### Individual d15N vs lamina layer, by species -----------------------------

## Least cisco

ggplot(data = lscs.delta.df, aes(x = laminae, y = d15N, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(4,14), breaks=seq(4,14,2), expand=c(0,0.5)) +
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, ncol=2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"),
    strip.text = element_text(size = 13),
    legend.position = "none",
    axis.line=element_line()
  )

## Humpback whitefish

ggplot(data = hbwf.delta.df, aes(x = laminae, y = d15N, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(3.5,12), breaks=seq(4,12,2), expand=c(0,0.5)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, ncol=2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"),
    strip.text = element_text(size = 13),
    legend.position = "none",
    axis.line=element_line()
  )

## Broad whitefish

ggplot(data = bdwf.delta.df, aes(x = laminae, y = d15N, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(3.5,12), breaks=seq(4,12,2), expand=c(0,0.5)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, ncol=2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"),
    strip.text = element_text(size = 13),
    legend.position = "none",
    axis.line=element_line()
  )

### Individual d13C vs lamina layer, by species -----------------------------

## Least cisco

ggplot(data = lscs.delta.df, aes(x = laminae, y = d13C, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_y_continuous(limits=c(-34,-20), breaks=seq(-34,-20, 2), expand=c(0,0.1)) +
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, ncol=2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"),
    strip.text = element_text(size = 13),
    legend.position = "none",
    axis.line=element_line()
  )

## Humpback whitefish

ggplot(data = hbwf.delta.df, aes(x = laminae, y = d13C, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_y_continuous(limits=c(-30,-19.5), breaks=seq(-30,-20, 2), expand=c(0,0.1)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, ncol=2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"),
    strip.text = element_text(size = 13),
    legend.position = "none",
    axis.line=element_line()
  )

## Broad whitefish

ggplot(data = bdwf.delta.df, aes(x = laminae, y = d13C, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_y_continuous(limits=c(-25,-20), breaks=seq(-24,-20, 2), expand=c(0,0.1)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, ncol=2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"),
    strip.text = element_text(size = 13),
    legend.position = "none",
    axis.line=element_line()
  )



## BDWF_45 two eye reconstruction ------------------------------------------
bdwf.45.df <- bdwf.delta.df %>% 
  filter(spp_ID %in% c("BDWF_45_1", "BDWF_45_2"))

ggplot(data = bdwf.45.df, aes(x = laminae, y = d13C, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_y_continuous(limits=c(-25,-20), breaks=seq(-25,-20, 1), minor_breaks = seq(-25,-20, 0.5), expand=c(0,0.1)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"),
    strip.text = element_text(size = 13),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.position = c(0.8,0.2),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.line=element_line()
  )

ggplot(data = bdwf.45.df, aes(x = laminae, y = d15N, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(4,11), breaks=seq(4,11,1), expand=c(0,0.1)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"),
    strip.text = element_text(size = 13),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.position = c(0.8,0.2),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.line=element_line()
  )


## Tissue isotope value correlations ---------------------------------------

### Muscle and lamina d15N correlations -------------------------------------

ggplot(correlation.df, aes(x = muscle.d15N, y = lamina.d15N, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  geom_point(shape = 17, size = 5) +
  xlab(expression(Muscle~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  ylab(expression(Lamina~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  scale_x_continuous(limits = c(8,14), breaks = seq(8,14,0.5), expand = c(0,0.1)) +
  scale_y_continuous(limits = c(8,14), breaks = seq(8,14,0.5), expand = c(0,0.1)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2")) +
  labs(color = "Sample Name") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    legend.position = c(0.1,0.7),
    legend.background = element_rect(fill = "gray90", color = "black"),
    legend.key = element_rect(fill = "gray90"),
    legend.title = element_text(size =14),
    legend.text = element_text(size=13),
    axis.line=element_line()
  )

muscle.lamina.d15N.lm <- lm(correlation.df$lamina.d15N ~ correlation.df$muscle.d15N)
summary(muscle.lamina.d15N.lm) # m = 0.45 (p = 0.01), b = 5.88 (p = <0.01), adj. R^2 = 0.89

# Same plot as above, just using the mean d15N for the final 3 lamina layers instead of just the 
# final layer:
#
# ggplot(correlation.df, aes(x = muscle.d15N, y = mean.lamina.d15N, color = Sample.Name)) +
#   geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
#   geom_point(shape = 17, size = 5) +
#   xlab(expression(Muscle~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
#   ylab(expression(Lamina~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
#   scale_x_continuous(limits = c(8,14), breaks = seq(8,14,0.5), expand = c(0,0.1)) +
#   scale_y_continuous(limits = c(8,14), breaks = seq(8,14,0.5), expand = c(0,0.1)) +
#   scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2")) +
#   labs(color = "Sample Name") +
#   theme(
#     panel.background = element_blank(),
#     axis.title.x = element_text(size=16, vjust = 0),
#     axis.title.y = element_text(size =16),
#     axis.text = element_text(size=13, color="black"), 
#     legend.position = c(0.1,0.7),
#     legend.background = element_rect(fill = "gray90", color = "black"),
#     legend.key = element_rect(fill = "gray90"),
#     legend.title = element_text(size =14),
#     legend.text = element_text(size=13),
#     axis.line=element_line()
#   )


### Muscle and lamina d13C correlations -------------------------------------

ggplot(correlation.df, aes(x = muscle.d13C, y = lamina.d13C, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  geom_point(shape = 16, size = 5) +
  xlab(expression(Muscle~italic(delta)^13*C~("\211"~VPDB))) +
  ylab(expression(Lamina~italic(delta)^13*C~("\211"~VPDB))) +
  scale_x_continuous(limits=c(-25,-19.5), breaks=seq(-25,-19.5, 0.5), expand=c(0,0.1)) +
  scale_y_continuous(limits=c(-25,-19.5), breaks=seq(-25,-19.5, 0.5), expand=c(0,0.1)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2")) +
  labs(color = "Sample Name") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    legend.position = "none",
    legend.background = element_rect(fill = "gray90", color = "black"),
    legend.key = element_rect(fill = "gray90"),
    legend.title = element_text(size =14),
    legend.text = element_text(size=13),
    axis.line=element_line()
  )

muscle.lamina.d13C.lm <- lm(correlation.df$lamina.d13C ~ correlation.df$muscle.d13C)
summary(muscle.lamina.d13C.lm) # m = -0.45 (p = 0.33), b = -33.31 (p = 0.03), adj. R^2 = 0.08

# Same plot as above, just using the mean d13C for the final 3 lamina layers instead of just the 
# final layer:
#
# ggplot(correlation.df, aes(x = muscle.d13C, y = mean.lamina.d13C, color = Sample.Name)) +
#   geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
#   geom_point(shape = 16, size = 5) +
#   xlab(expression(Muscle~italic(delta)^13*C~("\211"~VPDB))) +
#   ylab(expression(Lamina~italic(delta)^13*C~("\211"~VPDB))) +
#   scale_x_continuous(limits=c(-25,-19.5), breaks=seq(-25,-19.5, 0.5), expand=c(0,0.1)) +
#   scale_y_continuous(limits=c(-25,-19.5), breaks=seq(-25,-19.5, 0.5), expand=c(0,0.1)) +
#   scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2")) +
#   labs(color = "Sample Name") +
#   theme(
#     panel.background = element_blank(),
#     axis.title.x = element_text(size=16, vjust = 0),
#     axis.title.y = element_text(size =16),
#     axis.text = element_text(size=13, color="black"), 
#     legend.position = "none",
#     legend.background = element_rect(fill = "gray90", color = "black"),
#     legend.key = element_rect(fill = "gray90"),
#     legend.title = element_text(size =14),
#     legend.text = element_text(size=13),
#     axis.line=element_line()
#   )


### Muscle and fin d15N correlations ----------------------------------------

ggplot(correlation.df, aes(x = muscle.d15N, y = fin.d15N, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  geom_point(shape = 17, size = 5) +
  xlab(expression(Muscle~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  ylab(expression(Fin~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  scale_x_continuous(limits = c(8,14), breaks = seq(8,14,0.5), expand = c(0,0.1)) +
  scale_y_continuous(limits = c(8,14), breaks = seq(8,14,0.5), expand = c(0,0.1)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2")) +
  labs(color = "Sample Name") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    legend.position = c(0.1,0.7),
    legend.background = element_rect(fill = "gray90", color = "black"),
    legend.key = element_rect(fill = "gray90"),
    legend.title = element_text(size =14),
    legend.text = element_text(size=13),
    axis.line=element_line()
  )

muscle.fin.d15N.lm <- lm(fin.d15N ~ muscle.d15N, data = correlation.df)
summary(muscle.fin.d15N.lm) # m = 1.14 (p = <0.01), b = -1.29 (p = 0.37), adj. R^2 = 0.96


### Muscle and fin d13C correlations ----------------------------------------

ggplot(correlation.df, aes(x = muscle.d13C, y = fin.d13C, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  geom_point(shape = 16, size = 5) +
  xlab(expression(Muscle~italic(delta)^13*C~("\211"~VPDB))) +
  ylab(expression(Fin~italic(delta)^13*C~("\211"~VPDB))) +
  scale_x_continuous(limits=c(-25,-19.5), breaks=seq(-25,-19.5, 0.5), expand=c(0,0.1)) +
  scale_y_continuous(limits=c(-25,-19.5), breaks=seq(-25,-19.5, 0.5), expand=c(0,0.1)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2")) +
  labs(color = "Sample Name") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    legend.position = c(0.1,0.7),
    legend.background = element_rect(fill = "gray90", color = "black"),
    legend.key = element_rect(fill = "gray90"),
    legend.title = element_text(size =14),
    legend.text = element_text(size=13),
    axis.line=element_line()
  )

muscle.fin.d13C.lm <- lm(fin.d13C ~ muscle.d13C, data = correlation.df)
summary(muscle.fin.d13C.lm) # m = 1.28 (p = <0.01), b = 8.16 (p = 0.18), adj. R^2 = 0.91

fin.muscle.adjustment.factor <- mean(correlation.df$fin.d13C - correlation.df$muscle.d13C)
fin.muscle.adjustment.factor # -1.572


# URSA figures ------------------------------------------------------------

# average.delta.laminae <- data %>% 
#   group_by(species,laminae) %>% 
#   summarise(
#     "mean.d15N" = mean(d15N, na.rm = T),
#     "se.mean.d15N" = sqrt(var(d15N, na.rm = T)),
#     "mean.d13C" = mean(d13C, na.rm = T),
#     "se.mean.d13C" = sqrt(var(d13C, na.rm = T))
#   )
# average.delta.laminae
# 
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# # Average d15N per lamina 
# ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d15N, color = species)) +
#   geom_line(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
#   geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
#   ylab(expression(Average~italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
#   xlab("Lamina layer") +
#   scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
#   scale_fill_manual(values=cbPalette[c(1,6,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
#   scale_color_manual(values=cbPalette[c(1,6,2)]) +
#   scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1), expand = c(0,0.2)) +
#   scale_y_continuous(limits=c(6,13),breaks=seq(6,13,1), expand = c(0,0.1)) +
#   theme(
#     panel.background = element_blank(),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size =20),
#     axis.text = element_text(size=18, color="black"), 
#     legend.position = c(0.80, 0.22),
#     legend.key = element_rect(fill = "white"),
#     legend.text = element_text(size = 18),
#     legend.title = element_text(size = 20, face = "bold"),
#     axis.line=element_line()
#   )
# 
# ggsave("figures/URSA_plot_mean.d15n.vs.lamina.layer.png", device = "png", dpi = 800, width = 8, height = 6, units = "in")
# 
# 
# # Average d13C per lamina 
# ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d13C)) +
#   geom_line(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
#   geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
#   ylab(expression(Average~italic(delta)^13*C~("\211"~" VPDB"))) +
#   xlab("Lamina layer") +
#   scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
#   scale_fill_manual(values=cbPalette[c(1,6,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
#   scale_color_manual(values=cbPalette[c(1,6,2)]) +
#   scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1), expand = c(0,0.2)) +
#   scale_y_continuous(limits=c(-27,-21),breaks=seq(-27,-21,1), expand = c(0,0.2)) +
#   theme(
#     panel.background = element_blank(),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size =20),
#     axis.text = element_text(size=18, color="black"), 
#     legend.position = c(0.80, 0.22),
#     legend.key = element_rect(fill = "white"),
#     legend.text = element_text(size = 18),
#     legend.title = element_text(size = 20, face = "bold"),
#     axis.line=element_line()
#   )
# 
# ggsave("figures/URSA_plot_mean.d13c.vs.lamina.layer.png", device = "png", dpi = 800, width = 8, height = 6, units = "in")
# 
# 
# # Average d15N vs Average d13C, by lamina layer
# ggplot(data = average.delta.laminae, aes(x = mean.d13C, y = mean.d15N)) +
#   geom_path(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
#   geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
#   scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
#   xlab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
#   ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
#   scale_fill_manual(values=cbPalette[c(1,6,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
#   scale_color_manual(values=cbPalette[c(1,6,2)]) +
#   scale_x_continuous(limits=c(-27,-21),breaks=seq(-27,-21,1), expand = c(0,0.1)) +
#   scale_y_continuous(limits=c(6,12.6),breaks=seq(6,12,1), expand = c(0,0.1)) +
#   theme(
#     panel.background = element_blank(),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size =20),
#     axis.text = element_text(size=18, color="black"), 
#     legend.position = c(0.275, 0.885),
#     legend.background = element_blank(),
#     legend.key = element_rect(fill = "white"),
#     legend.text = element_text(size = 18),
#     legend.title = element_text(size = 20, face = "bold"),
#     axis.line=element_line()
#   )
# 
# ggsave("figures/URSA_plot_mean.d13c.vs.mean.d15n.png", device = "png", dpi = 800, width = 6, height = 6, units = "in")
# 
# 
# # Number of laminae by length
# length.dat <- read.csv("data/data.csv")
# 
# length.dat <- length.dat %>% 
#   unite("spp_ID", species:ID, sep= "_", remove = FALSE) %>% 
#   inner_join(data, length.dat, by = "spp_ID")
# 
# lin.model.dat <- length.dat %>% 
#   group_by(spp_ID) %>% 
#   summarize("species" = unique(species.x), "n.laminae" = length(laminae), "length_mm" = unique(length_mm))
# lin.model.dat
# lin.model.dat[3,3] <- 7
# lin.model.dat[4,3] <- 9
# lin.model <- lm(n.laminae ~ length_mm, data = lin.model.dat)
# summary(lin.model) # Slope significant (p-value = 0.004723); Adjusted R^2 = 0.3306
# 
# ggplot(data = lin.model.dat, aes(x = length_mm, y = n.laminae)) +
#   geom_abline(color = "black", cex = 1.2, slope = 0.028568, intercept = -1.875996) +
#   geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
#   scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
#   scale_fill_manual(values=cbPalette[c(1,6,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
#   scale_x_continuous(limits=c(220,400),breaks=seq(225,400,25), expand = c(0,5)) +
#   scale_y_continuous(limits=c(2.8,11.2),breaks=seq(3,11,1), expand = c(0,0.1)) +
#   xlab("Fish length (mm)") +
#   ylab("Number of laminae") +
#   theme(
#     panel.background = element_blank(),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size =20),
#     axis.text = element_text(size=18, color="black"), 
#     legend.position = c(0.22, 0.85),
#     legend.key = element_rect(fill = "white"),
#     legend.text = element_text(size = 18),
#     legend.title = element_text(size = 20, face = "bold"),
#     axis.line=element_line()
#   )
# 
# ggsave("figures/URSA_plot_length.vs.n.lamina.png", device = "png", dpi = "retina", width = 6, height = 6, units = "in")
# 
# 
# 
# final.df <- data %>% 
#   filter(spp_ID %in% c("BDWF_45_1", "HBWF_39", "LSCS_38"))
# final.df  
# 
# id.labs <- c("Broad Whitefish ID #45","Humpback Whitefish ID #39","Least Cisco ID #39")
# names(id.labs) <- c("BDWF_45_1", "HBWF_39", "LSCS_38")
# 
# ggplot(data = final.df, aes(x = laminae, y = d15N)) +
#   geom_line(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
#   geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
#   scale_shape_manual(values=c(21,22,24))+
#   scale_fill_manual(values=cbPalette[c(1,6,2)]) +
#   scale_color_manual(values=cbPalette[c(1,6,2)]) +
#   xlab("Lamina layer") +
#   ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
#   scale_y_continuous(limits=c(4,13), breaks=seq(4,13,1), expand=c(0,0.2)) +
#   scale_x_continuous(limits = c(0,11), breaks = seq(0,11,1), expand = c(0,0.5)) +
#   facet_wrap(~spp_ID, nrow = 1, scales="free_y", labeller = labeller(spp_ID = id.labs)) +
#   theme(
#     panel.background = element_blank(),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size =20),
#     axis.text = element_text(size=18, color="black"), 
#     strip.text = element_text(size = 18),
#     legend.position = "none",
#     axis.line=element_line()
#   )
# 
# ggsave("figures/example.d15n.vs.lamina.layer.png", device = "png", dpi = 800, width = 11.5, height = 4.2, units = "in")



# C:N ratios --------------------------------------------------------------


# Carbon 2 end-member mixing model ----------------------------------------


