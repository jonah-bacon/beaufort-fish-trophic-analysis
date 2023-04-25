# =====================================================================
# Lens chronology analysis
# Jonah Bacon
# 27 February 2023
# =====================================================================


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggpmisc)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

species.names <- c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")
names(species.names) <- c("ARCS","BDWF","HBWF","LSCS")


# Load data ---------------------------------------------------------------

lens.chrono <- read.csv("data/data - Lens Chronology - Sectioned.csv")
lens.chrono <- lens.chrono %>% 
  select(spp_ID, species, ID, Measurement_um, N_layers_LD)

lw.df <- read.csv("data/data - Sample Data.csv")
lw.df <- lw.df %>% select(spp_ID, species, ID, length_mm, weight_kg)

otolith.df <- read.csv("data/data - Otolith Chronology.csv", header = TRUE)
otolith.df <- otolith.df %>% 
  select(spp_ID, species, ID, age_estimate_LD, age_estimate_SB) %>% 
  group_by(spp_ID, species, ID) %>% 
  summarize(mean.age = ifelse(is.na(age_estimate_LD), age_estimate_SB, mean(age_estimate_LD, age_estimate_SB)))

lens.chrono.df <- right_join(lw.df, lens.chrono, by = c("spp_ID", "species", "ID"))
complete.df <- left_join(lens.chrono.df, otolith.df, by = c("spp_ID", "species", "ID"))

# Visualize data ----------------------------------------------------------

ggplot(complete.df, aes(x = mean.age, y = length_mm)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab("Fork length (mm)") +
  xlab("Age") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_x_continuous(limits=c(-0.5,15.5),breaks=seq(0,15,3), expand = c(0,0)) +
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
# ggsave("figures/length-at-age.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(complete.df, aes(x = mean.age, y = Measurement_um, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab("Lens diameter (\u03bcm)") +
  xlab("Age") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_x_continuous(limits=c(-0.5,15.5),breaks=seq(0,15,3), expand = c(0,0)) +
  scale_y_continuous(limits=c(450,3050),breaks=seq(500,3000,500), expand = c(0,0)) +
  facet_wrap(~species, nrow = 2, labeller = labeller(species = species.names)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label))), color = "black", label.x = 0.95, label.y = 0.95, parse = TRUE) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(rr.label))), color = "black", label.x = 0.95, label.y = 0.85, parse = TRUE) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=10, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave("figures/lens.diam.vs.age.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(data = complete.df, aes(x = length_mm, y = Measurement_um, color = species)) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  xlab("Fork length (mm)") +
  ylab("Lens diameter (\u03bcm)") +
  scale_shape_manual(values=c(21,22,24,25), name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)]) +
  scale_x_continuous(limits=c(80,540),breaks=seq(100,500,100), expand = c(0,0)) +
  scale_y_continuous(limits=c(450,3050),breaks=seq(500,3000,500), expand = c(0,0)) +
  facet_wrap(~species, nrow = 2, labeller = labeller(species = species.names)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label))), color = "black", label.x = 0.05, label.y = 0.95, parse = TRUE) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(rr.label))), color = "black", label.x = 0.05, label.y = 0.85, parse = TRUE) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )
# ggsave("figures/lens.diam.vs.FL.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(data = filter(complete.df, !is.na(N_layers_LD)), aes(y = mean.age, x = as.factor(N_layers_LD), color = species)) +
  geom_boxplot(aes(fill = species), color = "black") +
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lty = 2, color = "gray40") +
  xlab("Number of lamina layers") +
  ylab("Age") +
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_x_continuous(limits=c(0.5,540),breaks=seq(100,500,100), expand = c(0,0)) +
  scale_y_continuous(limits=c(-0.5,15.5),breaks=seq(0,15,3), expand = c(0,0)) +
  facet_wrap(~species, nrow = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = "none",
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    axis.line=element_line()
  )
# ggsave("figures/age.vs.nlamina.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(data = filter(complete.df, !is.na(N_layers_LD)), aes(y = length_mm, x = as.factor(N_layers_LD), color = species)) +
  geom_boxplot(aes(fill = species), color = "black") +
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lty = 2, color = "gray40") +
  xlab("Number of lamina layers") +
  ylab("Fork length (mm)") +
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_x_continuous(limits=c(0.5,540),breaks=seq(100,500,100), expand = c(0,0)) +
  scale_y_continuous(limits=c(80,520),breaks=seq(100,500,100), expand = c(0,0)) +
  facet_wrap(~species, nrow = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = "none",
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    axis.line=element_line()
  )
# ggsave("figures/FL.vs.nlamina.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

ggplot(data = filter(complete.df, !is.na(N_layers_LD)), aes(y = Measurement_um, x = as.factor(N_layers_LD))) +
  geom_boxplot(aes(fill = species), color = "black") +
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lty = 2, color = "gray40") +
  xlab("Number of lamina layers") +
  ylab("Lens diameter (\u03bcm)") +
  scale_fill_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2,4)], name = "Species", labels = c("Arctic Cisco", "Broad Whitefish", "Humpback Whitefish", "Least Cisco")) +
  # scale_x_continuous(limits=c(0.5,540),breaks=seq(100,500,100), expand = c(0,0)) +
  # scale_y_continuous(limits=c(80,520),breaks=seq(100,500,100), expand = c(0,0)) +
  facet_wrap(~species, nrow = 2, labeller = labeller(species = species.names)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.key = element_rect(fill = NA),
    legend.position = "none",
    legend.title = element_text(face = "bold.italic", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "gray80", color = "black"),
    axis.line=element_line()
  )
# ggsave("figures/lens.diam.vs.nlamina.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")


# Isometric lens, length, age relationships -------------------------------

arcs.lm <- lm(length_mm ~ Measurement_um, complete.df, subset = species == "ARCS")
bdwf.lm <- lm(length_mm ~ Measurement_um, complete.df, subset = species == "BDWF")
hbwf.lm <- lm(length_mm ~ Measurement_um, complete.df, subset = species == "HBWF")
lscs.lm <- lm(length_mm ~ Measurement_um, complete.df, subset = species == "LSCS")

load("data/bulk_sia_data.RData")

temp.df <- data.frame(ARCS = rep(NA, 119), BDWF = rep(NA, 119), HBWF = rep(NA, 119), LSCS = rep(NA, 119))
i=1

  df <- data.frame(filter(sia.df, species == names(species.names)[i]))
arcs.preds <- data.frame(SIAposition = filter(sia.df, species == "ARCS")$SIAposition, iso.length = predict(arcs.lm, newdata = data.frame(Measurement_um = filter(sia.df, species == "ARCS")$SIAposition)))
bdwf.preds <- predict(bdwf.lm, newdata = data.frame(Measurement_um = filter(sia.df, species == "BDWF")$SIAposition))
hbwf.preds <- predict(hbwf.lm, newdata = data.frame(Measurement_um = filter(sia.df, species == "HBWF")$SIAposition))
lscs.preds <- predict(lscs.lm, newdata = data.frame(Measurement_um = filter(sia.df, species == "LSCS")$SIAposition))


sia.df %>% 
  group_by(species) %>% 
  summarise(N = n())
  summarize("iso.length" = predict.lm(arcs.lm, newdata = data.frame(Measurement_um = df$SIAposition)))
predict(arcs.lm, class(complete.df$Measurement_um))
