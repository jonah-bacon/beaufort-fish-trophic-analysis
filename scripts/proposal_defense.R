# =====================================================================
# Bulk carbon & nitrogen SIA from eye lenses
# Jonah Bacon
# 25 April 2022
# =====================================================================

# Load packages -----------------------------------------------------------

library(ggplot2)
library(tidyverse)
# library(ggsci)
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
  select(-c("V9", "V10", "V11", "V12", "V13", "V14")) 
colnames(muscle.fin.df) <- c("Sample.Name","Sample.wt.mg","N.sig.V","C.sig.V","Conc.N","Conc.C","d15N","d13C")

muscle.fin.df <- muscle.fin.df %>% 
  mutate("C:N" = Conc.C/Conc.N)

standards.df <- muscle.fin.df %>% 
  filter(Sample.Name == "ref/chk/peptone")

muscle.fin.df <- muscle.fin.df %>% 
  filter(Sample.Name != "ref/chk/peptone") %>% 
  separate("Sample.Name", into = c("Sample.Name", "Sample.Type"), sep = 4, remove = T)

muscle.fin.df$Sample.Name <- c("BDWF_41", "LSCS_42", "HBWF_33", "LSCS_44", "BDWF_45_1", "LSCS_42", "HBWF_33", "LSCS_44", "BDWF_41", "BDWF_45_1")

muscle.fin.df$Sample.Name <- as.factor(muscle.fin.df$Sample.Name)
muscle.fin.df$Sample.Type <- as.factor(muscle.fin.df$Sample.Type)

str(muscle.fin.df)

muscle.fin.df1 <- muscle.fin.df %>% 
  filter(Sample.Type == "F") %>% 
  mutate(Sample.Type = "Fin") %>% 
  summarise("Sample.Name" = Sample.Name, "fin.d15N" = d15N, "fin.d13C" = d13C)

muscle.fin.df2 <- muscle.fin.df %>% 
  filter(Sample.Type == "M") %>% 
  mutate(Sample.Type = "Muscle") %>% 
  summarise("Sample.Name" = Sample.Name, "muscle.d15N" = d15N, "muscle.d13C" = d13C)

# muscle.fin.df <- rbind(muscle.fin.df1, muscle.fin.df2)
# muscle.fin.df <- muscle.fin.df %>% 
#   select(Sample.Name, Sample.Type, d15N, d13C)
# muscle.fin.df

temp.df <- delta.df %>% 
  select(spp_ID, laminae, d15N, d13C) %>% 
  filter(spp_ID %in% c("BDWF_41", "LSCS_42", "HBWF_33", "LSCS_44", "BDWF_45_1")) %>% 
  mutate(Sample.Name = spp_ID, Sample.Type = laminae, d15N = d15N, d13C = d13C, .keep = "none") %>% 
  group_by(Sample.Name) %>% 
  summarise("lamina.d15N" = d15N[1], "lamina.d13C" = d13C[1], "mean.lamina.d15N" = mean(d15N[1:3]), "mean.lamina.d13C" = mean(d13C[1:3]))


combined.df <- merge(muscle.fin.df1, muscle.fin.df2)
combined.df <- merge(combined.df, temp.df)

# Visualize the data ------------------------------------------------------

## d13C vs lamina layer, faceted by species

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

## d15N vs lamina layer, faceted by species

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

# Fin vs muscle vs lens SIA -----------------------------------------------

ggplot(combined.df, aes(x = muscle.d15N, y = lamina.d15N, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  # geom_abline(slope = 1, intercept = d15N.intercept.hat, cex = 1.5, color = "black") +
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

ggplot(combined.df, aes(x = muscle.d15N, y = mean.lamina.d15N, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  # geom_abline(slope = 1, intercept = d15N.intercept.hat, cex = 1.5, color = "black") +
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

ggplot(combined.df, aes(x = muscle.d13C, y = lamina.d13C, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  # geom_abline(slope = 1, intercept = d15N.intercept.hat, cex = 1.5, color = "black") +
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

ggplot(combined.df, aes(x = muscle.d13C, y = lamina.d13C, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  # geom_abline(slope = 1, intercept = d13C.intercept.hat, color = "black", cex = 1.5) +
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

ggplot(combined.df, aes(x = muscle.d13C, y = mean.lamina.d13C, color = Sample.Name)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
  # geom_abline(slope = 1, intercept = d13C.intercept.hat, color = "black", cex = 1.5) +
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
