### Bulk C/N Stable Isotope Laminae Analysis

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)

rough.df <- read.csv("data/first_bulk_data.csv", header = T)

head(rough.df)
str(rough.df)

not.so.rough.df <- rough.df %>% select(-c(layer, data_sheet_ID, well_ID, tray_ID, Sample.Name, Sample.Wt., N.Signal, C.Signal, Conc.C, Conc.N))

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
