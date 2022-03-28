### Bulk C/N Stable Isotope Laminae Analysis

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)
library(ggrepel)

rough.df <- read.csv("data/data - Bulk SIA Data.csv", header = T)

head(rough.df)
str(rough.df)

not.so.rough.df <- rough.df %>% select(-c(layer, data_sheet_ID, well_ID, tray_ID, Sample.Name, Sample.Wt., N.Signal, C.Signal, Conc.C, Conc.N, Notes))

data <- not.so.rough.df %>%
  unite("spp_ID", species:ID, sep= "_", remove = FALSE)

data$spp_ID <- as.factor(data$spp_ID)
data$species <- as.factor(data$species)
data$ID <- as.factor(data$ID)
        
head(data)
str(data)

lscs.df <- data %>% filter(species == "LSCS")
bdwf.df <- data %>% filter(species == "BDWF")
hbwf.df <- data %>% filter(species == "HBWF")


# Visualize d15N data

ggplot(data = data, aes(x = laminae, y = d15N), ylab("Laminae layer")) +
  geom_point(aes(col = species)) +
  xlab("Laminae layer")

ggplot(data = data, aes(x = laminae, y = d15N, col = spp_ID)) +
  geom_line() +
  geom_point() +
  xlab("Laminae layer") +
  labs(col = "Spp. + ID")

### Least Cisco
ggplot(data = lscs.df, aes(x = laminae, y = d15N, col = ID)) +
  geom_line() +
  geom_point() +
  xlab("Laminae layer") +
  labs(col = "ID")

### Humpback Whitefish
ggplot(data = hbwf.df, aes(x = laminae, y = d15N, col = ID)) +
  geom_line() +
  geom_point() +
  xlab("Laminae layer") +
  labs(col = "ID")

### Broad Whitefish
ggplot(data = bdwf.df, aes(x = laminae, y = d15N, col = ID)) +
  geom_line() +
  geom_point() +
  xlab("Laminae layer") +
  labs(col = "ID")

ggplot(data = data, aes(x=laminae, y=d15N, color=spp_ID)) +
  geom_line() +
  geom_point() +
  facet_wrap(~spp_ID, ncol=2)

i=1
for (i in 1:length(unique(data$spp_ID))) {
  test.plot <- ggplot(data = data[data$spp_ID == levels(data$spp_ID)[i],], aes(x=laminae, y=d15N)) +
    geom_line() +
    geom_point() +
    labs(title = levels(data$spp_ID)[i])
  print(test.plot)
}
### Now you can click through the plots in that lower right-hand plot window - There should be 11 of them, 1 for each individual

# Visualize d13C data

ggplot(data = data, aes(x = laminae, y = d13C), ylab("Laminae layer")) +
  geom_point(aes(col = species)) +
  xlab("Laminae layer")

ggplot(data = data, aes(x = laminae, y = d13C, col = spp_ID)) +
  geom_line() +
  geom_point() +
  xlab("Laminae layer") +
  labs(col = "Spp. + ID")

### Least Cisco
ggplot(data = lscs.df, aes(x = laminae, y = d13C, col = ID)) +
  geom_line() +
  geom_point() +
  xlab("Laminae layer") +
  labs(col = "ID")

### Humpback Whitefish
ggplot(data = hbwf.df, aes(x = laminae, y = d13C, col = ID)) +
  geom_line() +
  geom_point() +
  xlab("Laminae layer") +
  labs(col = "ID")

### Broad Whitefish
ggplot(data = bdwf.df, aes(x = laminae, y = d13C, col = ID)) +
  geom_line() +
  geom_point() +
  xlab("Laminae layer") +
  labs(col = "ID")

ggplot(data = data, aes(x=laminae, y=d13C, color=spp_ID)) +
  geom_line() +
  geom_point() +
  facet_wrap(~spp_ID, ncol=2)

i=1
for (i in 1:11) {
  test.plot <- ggplot(data = data[data$spp_ID == levels(data$spp_ID)[i],], aes(x=laminae, y=d13C)) +
    geom_line() +
    geom_point() +
    labs(title = levels(data$spp_ID)[i])
  print(test.plot)
}


# AFS figures -------------------------------------------------------------

afs.df <- hbwf.df %>% filter(ID %in% c("31","36","38","39"))

id.labs <- c("Humpback Whitefish ID #31","Humpback Whitefish ID #36","Humpback Whitefish ID #38","Humpback Whitefish ID #39")
names(id.labs) <- c("HBWF_31","HBWF_36","HBWF_38","HBWF_39")

ggplot(data = afs.df, aes(x = laminae, y = d15N, color = spp_ID)) +
  geom_point(cex = 5) +
  geom_line(cex = 2, alpha = 0.5) +
  scale_fill_jco() +
  scale_color_jco() +
  geom_point(shape = 21, color = "black", cex = 5) +
  xlab("Laminae Layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(3.5,12.5), breaks=seq(4,12,1), expand=c(0,0)) +
  scale_x_continuous(limits = c(0,11), breaks = seq(0,11,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, ncol=2, scales="free_y", labeller = labeller(spp_ID = id.labs)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    strip.text = element_text(size = 13),
    legend.position = "none",
    axis.line=element_line()
  )

ggsave("HBWF_AFS_plot.png", device = "png", dpi = 800)




# URSA Figures ------------------------------------------------------------


average.delta.laminae <- data %>% 
  group_by(species,laminae) %>% 
  summarise(
    "mean.d15N" = mean(d15N, na.rm = T),
    "se.mean.d15N" = sqrt(var(d15N, na.rm = T)),
    "mean.d13C" = mean(d13C, na.rm = T),
    "se.mean.d13C" = sqrt(var(d13C, na.rm = T))
  )
average.delta.laminae

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Average d15N per lamina 
ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d15N, color = species)) +
  geom_line(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab(expression(Average~italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  xlab("Lamina layer") +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  scale_fill_manual(values=cbPalette[c(1,6,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,6,2)]) +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1), expand = c(0,0.2)) +
  scale_y_continuous(limits=c(6,13),breaks=seq(6,13,1), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size =20),
    axis.text = element_text(size=18, color="black"), 
    legend.position = c(0.80, 0.22),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.line=element_line()
  )

ggsave("figures/URSA_plot_mean.d15n.vs.lamina.layer.png", device = "png", dpi = 800, width = 8, height = 6, units = "in")


# Average d13C per lamina 
ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d13C)) +
  geom_line(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab(expression(Average~italic(delta)^13*C~("\211"~" VPDB"))) +
  xlab("Lamina layer") +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  scale_fill_manual(values=cbPalette[c(1,6,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,6,2)]) +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1), expand = c(0,0.2)) +
  scale_y_continuous(limits=c(-27,-21),breaks=seq(-27,-21,1), expand = c(0,0.2)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size =20),
    axis.text = element_text(size=18, color="black"), 
    legend.position = c(0.80, 0.22),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.line=element_line()
  )

ggsave("figures/URSA_plot_mean.d13c.vs.lamina.layer.png", device = "png", dpi = 800, width = 8, height = 6, units = "in")


# Average d15N vs Average d13C, by lamina layer
ggplot(data = average.delta.laminae, aes(x = mean.d13C, y = mean.d15N)) +
  geom_path(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  xlab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_fill_manual(values=cbPalette[c(1,6,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,6,2)]) +
  scale_x_continuous(limits=c(-27,-21),breaks=seq(-27,-21,1), expand = c(0,0.1)) +
  scale_y_continuous(limits=c(6,12.6),breaks=seq(6,12,1), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size =20),
    axis.text = element_text(size=18, color="black"), 
    legend.position = c(0.275, 0.885),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.line=element_line()
)

ggsave("figures/URSA_plot_mean.d13c.vs.mean.d15n.png", device = "png", dpi = 800, width = 6, height = 6, units = "in")


# Number of laminae by length
length.dat <- read.csv("data/data.csv")

length.dat <- length.dat %>% 
  unite("spp_ID", species:ID, sep= "_", remove = FALSE) %>% 
  inner_join(data, length.dat, by = "spp_ID")

lin.model.dat <- length.dat %>% 
  group_by(spp_ID) %>% 
  summarize("species" = unique(species.x), "n.laminae" = length(laminae), "length_mm" = unique(length_mm))
lin.model.dat
lin.model.dat[3,3] <- 7
lin.model.dat[4,3] <- 9
lin.model <- lm(n.laminae ~ length_mm, data = lin.model.dat)
summary(lin.model) # Slope significant (p-value = 0.004723); Adjusted R^2 = 0.3306

ggplot(data = lin.model.dat, aes(x = length_mm, y = n.laminae)) +
  geom_abline(color = "black", cex = 1.2, slope = 0.028568, intercept = -1.875996) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,6,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_x_continuous(limits=c(220,400),breaks=seq(225,400,25), expand = c(0,5)) +
  scale_y_continuous(limits=c(2.8,11.2),breaks=seq(3,11,1), expand = c(0,0.1)) +
  xlab("Fish length (mm)") +
  ylab("Number of laminae") +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size =20),
    axis.text = element_text(size=18, color="black"), 
    legend.position = c(0.22, 0.85),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.line=element_line()
  )

ggsave("figures/URSA_plot_length.vs.n.lamina.png", device = "png", dpi = "retina", width = 6, height = 6, units = "in")



final.df <- data %>% 
  filter(spp_ID %in% c("BDWF_45_1", "HBWF_39", "LSCS_38"))
final.df  

id.labs <- c("Broad Whitefish ID #45","Humpback Whitefish ID #39","Least Cisco ID #39")
names(id.labs) <- c("BDWF_45_1", "HBWF_39", "LSCS_38")

ggplot(data = final.df, aes(x = laminae, y = d15N)) +
  geom_line(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  scale_shape_manual(values=c(21,22,24))+
  scale_fill_manual(values=cbPalette[c(1,6,2)]) +
  scale_color_manual(values=cbPalette[c(1,6,2)]) +
  xlab("Lamina layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(4,13), breaks=seq(4,13,1), expand=c(0,0.2)) +
  scale_x_continuous(limits = c(0,11), breaks = seq(0,11,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, nrow = 1, scales="free_y", labeller = labeller(spp_ID = id.labs)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size =20),
    axis.text = element_text(size=18, color="black"), 
    strip.text = element_text(size = 18),
    legend.position = "none",
    axis.line=element_line()
  )

ggsave("figures/example.d15n.vs.lamina.layer.png", device = "png", dpi = 800, width = 11.5, height = 4.2, units = "in")
