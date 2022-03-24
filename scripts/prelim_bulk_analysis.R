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
for (i in 1:11) {
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





average.delta.laminae <- data %>% 
  group_by(species,laminae) %>% 
  summarise(
    "mean.d15N" = mean(d15N, na.rm = T),
    "se.mean.d15N" = sqrt(var(d15N, na.rm = T)),
    "mean.d13C" = mean(d13C, na.rm = T),
    "se.mean.d13C" = sqrt(var(d13C, na.rm = T))
  )
average.delta.laminae

ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d15N, color = species)) +
  geom_point(cex = 3, position=position_dodge(0.2)) +
  geom_line(cex = 0.8) +
  geom_errorbar(aes(ymin=mean.d15N-se.mean.d15N, ymax=mean.d15N+se.mean.d15N), cex = 0.8, width=.2, position=position_dodge(0.2)) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  xlab("Lamina layer") +
  labs(color = "Species") +
  scale_x_continuous(limits=c(-0.2,10.2),breaks=seq(0,11,1), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    strip.text = element_text(size = 13),
    legend.position = c(0.7, 0.2),
    axis.line=element_line()
  )

# Average d15N per lamina 
ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d15N, color = species)) +
  geom_ribbon(aes(ymin = mean.d15N-se.mean.d15N, ymax = mean.d15N+se.mean.d15N, fill = species), alpha = 0.2, show.legend = F) +
  geom_line(cex = 1.2) +
  geom_point(cex = 4) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  xlab("Lamina layer") +
  scale_color_discrete(name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    legend.position = c(0.8, 0.2),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.line=element_line()
  )

# Average d13C per lamina 
ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d13C, color = species)) +
  geom_ribbon(aes(ymin = mean.d13C-se.mean.d13C, ymax = mean.d13C+se.mean.d13C, fill = species), alpha = 0.2, show.legend = F) +
  geom_line(cex = 1.2) +
  geom_point(cex = 4) +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  xlab("Lamina layer") +
  scale_color_discrete(name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,11,1), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    legend.position = c(0.8, 0.2),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.line=element_line()
  )

average.delta.laminae

# Combined plot
ggplot(data = average.delta.laminae, aes(x = mean.d13C, y = mean.d15N, color = species, label = laminae)) +
  geom_path(cex = 1.2, alpha = 0.5) +
  geom_errorbar(aes(ymin=mean.d15N-se.mean.d15N, ymax=mean.d15N+se.mean.d15N), cex = 0.5, width=.1) +
  geom_errorbar(aes(xmin = mean.d13C-se.mean.d13C, xmax = mean.d13C+se.mean.d13C), cex = 0.5, width=.1, key_glyph = "point") +
  geom_point(aes(fill = species), shape = 21, color = "black", cex = 5, show.legend = F) +
  # geom_label_repel(box.padding = 1, color = "black") +
  xlab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_color_discrete(name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_x_continuous(limits=c(-29.5,-20.5),breaks=seq(-29,-21,1), expand = c(0,0.1)) +
  scale_y_continuous(limits=c(3.5,12.5),breaks=seq(4,12,1), expand = c(0,0.1)) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))+
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=16, vjust = 0),
    axis.title.y = element_text(size =16),
    axis.text = element_text(size=13, color="black"), 
    legend.position = c(0.18, 0.14),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.line=element_line()
  )

?draw_key
