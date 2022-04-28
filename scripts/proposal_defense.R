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

# {240,163,255},{0,117,220},{153,63,0},{76,0,92},{25,25,25},{0,92,49},{43,206,72},{255,204,153},{128,128,128},
# {148,255,181},{143,124,0},{157,204,0},{194,0,136},{0,51,128},{255,164,5},{255,168,187},{66,102,0},{255,0,16},
# {94,241,242},{0,153,143},{224,255,102},{116,10,255},{153,0,0},{255,255,128},{255,255,0},{255,80,5}.

color_pal <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
"#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
"#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
"#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
"#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
"#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
"#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
"#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
"#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
"#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
"#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
"#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
"#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
"#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
"#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
"#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
"#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
"#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
"#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
"#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
"#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
"#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"
)

color_pal1 <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00")
color_pal2 <- c("#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
 # orange
  "black",  "gold",
  "gray50", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1",
  "darkorange4", "brown"
)

delta.df %>% 
  group_by(ID) %>% 
  summarise("n" = n())
#f
# Load eye lens data ------------------------------------------------------

bulk.data <- read.csv("data/bulk.SIA.data.csv", header = T)

head(bulk.data)
tail(bulk.data)
str(bulk.data)

delta.df <- bulk.data %>% 
  select(c(species, ID, laminae, d15N, d13C)) %>% 
  unite("spp_ID", species:ID, sep= "_", remove = FALSE) %>% 
  na.omit()

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

average.delta.laminae <- delta.df %>%
  group_by(species,laminae) %>%
  summarise(
    "mean.d15N" = mean(d15N, na.rm = T),
    "se.mean.d15N" = sqrt(var(d15N, na.rm = T)),
    "mean.d13C" = mean(d13C, na.rm = T),
    "se.mean.d13C" = sqrt(var(d13C, na.rm = T))
  )
average.delta.laminae

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

correlation.df <- merge(
                     delta.df %>% 
                       select(species, spp_ID, laminae, d15N, d13C) %>% 
                       filter(spp_ID %in% c("BDWF_41", "LSCS_42", "HBWF_33", "LSCS_44", "BDWF_45_1")) %>% 
                       mutate(species = species, Sample.Name = spp_ID, Sample.Type = laminae, d15N = d15N, d13C = d13C, .keep = "none") %>% 
                       group_by(species, Sample.Name) %>% 
                       summarise("lamina.d15N" = d15N[1], "lamina.d13C" = d13C[1], "mean.lamina.d15N" = mean(d15N[1:3]), "mean.lamina.d13C" = mean(d13C[1:3])),
                     muscle.fin.df
                     )

# Visualize the data ------------------------------------------------------

## d13C vs lamina layer, faceted by species --------------------------------

ggplot(data = delta.df, aes(x = laminae, y = d13C)) +
  geom_line(aes(col = ID), cex = 1.2, alpha = 0.5) +
  geom_point(aes(col = ID), cex = 3) +
  xlab("Lamina layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_color_simpsons() +
  scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1), expand=c(0,0.1)) +
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

# ggsave("figures/species-faceted.d13C.per.lamina.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

## d15N vs lamina layer, faceted by species --------------------------------

ggplot(data = delta.df, aes(x = laminae, y = d15N)) +
  geom_line(aes(col = ID), cex = 1.2, alpha = 0.7) +
  geom_point(aes(col = ID), cex = 3) +
  xlab("Lamina layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_color_simpsons() +
  scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1), expand=c(0,0.1)) +
  scale_y_continuous(limits = c(3.5,13), breaks = seq(4,12,2), expand = c(0,0.1)) +
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

# ggsave("figures/species-faceted.d15N.per.lamina.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

# d15N.df <- fish.df %>%
#   select(Sample.Name, Sample.Type, d15N) %>% 
#   pivot_wider(names_from = Sample.Type, values_from = d15N)
# 
# d15N.intercept.hat <- mean(d15N.df$F - d15N.df$M)
# d15N.intercept.hat
# 
# d15N.lm <- lm(F ~ M, data = d15N.df)
# summary(d15N.lm)
# 
# ggplot(d15N.df, aes(x = Muscle, y = Fin, color = Sample.Name)) +
#   geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, alpha = 0.4) +
#   geom_abline(slope = 1, intercept = d15N.intercept.hat, cex = 1.5, color = "black") +
#   geom_point(shape = 17, size = 5) +
#   xlab(expression(Muscle~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
#   ylab(expression(Fin~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
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


## d15N vs lamina layer, by species and ID -----------------------------

## Least cisco

ggplot(data = lscs.delta.df, aes(x = laminae, y = d15N, color = spp_ID)) +
  geom_line(cex = 0.8, color = "gray80") +
  geom_point(cex = 2, color = "black") +
  xlab("Lamina Layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(4,14), breaks=seq(4,14,2), expand=c(0,0.5)) +
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, nrow = 2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=12, vjust = 0),
    axis.title.y = element_text(size =12),
    axis.text = element_text(size=8, color="black"),
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/individual.d15N.per.lamina.LSCS.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Humpback whitefish

ggplot(data = hbwf.delta.df, aes(x = laminae, y = d15N, color = spp_ID)) +
  geom_line(cex = 0.8, color = "gray80") +
  geom_point(cex = 2, color = "black") +
  xlab("Lamina Layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(3.5,12), breaks=seq(4,12,2), expand=c(0,0.5)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, nrow = 2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=12, vjust = 0),
    axis.title.y = element_text(size =12),
    axis.text = element_text(size=8, color="black"),
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/individual.d15N.per.lamina.HBWF.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Broad whitefish

ggplot(data = bdwf.delta.df, aes(x = laminae, y = d15N, color = spp_ID)) +
  geom_line(cex = 0.8, color = "gray80") +
  geom_point(cex = 2, color = "black") +
  xlab("Lamina Layer") +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]~" gas"))) +
  scale_y_continuous(limits=c(3.5,12), breaks=seq(4,12,2), expand=c(0,0.5)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, nrow = 2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=12, vjust = 0),
    axis.title.y = element_text(size =12),
    axis.text = element_text(size=8, color="black"),
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/individual.d15N.per.lamina.BDWF.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## d13C vs lamina layer, by species and ID -----------------------------

## Least cisco

ggplot(data = lscs.delta.df, aes(x = laminae, y = d13C, color = spp_ID)) +
  geom_line(cex = 0.8, color = "gray80") +
  geom_point(cex = 2, color = "black") +
  xlab("Lamina Layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_y_continuous(limits=c(-35,-20), breaks=seq(-35,-20, 3), expand=c(0,0.1)) +
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, nrow = 2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=12, vjust = 0),
    axis.title.y = element_text(size =12),
    axis.text = element_text(size=8, color="black"),
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/individual.d13C.per.lamina.LSCS.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Humpback whitefish

ggplot(data = hbwf.delta.df, aes(x = laminae, y = d13C, color = spp_ID)) +
  geom_line(cex = 0.8, color = "gray80") +
  geom_point(cex = 2, color = "black") +
  xlab("Lamina Layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_y_continuous(limits=c(-30,-19.5), breaks=seq(-30,-20, 2), expand=c(0,0.1)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, nrow = 2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=12, vjust = 0),
    axis.title.y = element_text(size =12),
    axis.text = element_text(size=8, color="black"),
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/individual.d13C.per.lamina.HBWF.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Broad whitefish

ggplot(data = bdwf.delta.df, aes(x = laminae, y = d13C, color = spp_ID)) +
  geom_line(cex = 0.8, color = "gray80") +
  geom_point(cex = 2, color = "black") +
  xlab("Lamina Layer") +
  ylab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  scale_y_continuous(limits=c(-25,-20), breaks=seq(-24,-20, 2), expand=c(0,0.1)) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), expand = c(0,0.5)) +
  facet_wrap(~spp_ID, nrow = 2, scales="free_y") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size=12, vjust = 0),
    axis.title.y = element_text(size =12),
    axis.text = element_text(size=8, color="black"),
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/individual.d13C.per.lamina.BDWF.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## BDWF_45 two eye reconstruction ------------------------------------------
bdwf.45.df <- bdwf.delta.df %>% 
  filter(spp_ID %in% c("BDWF_45_1", "BDWF_45_2"))

ggplot(data = bdwf.45.df, aes(x = laminae, y = d13C, color = spp_ID, fill = spp_ID)) +
  geom_line(cex = 1.2, alpha = 0.5) +
  geom_point(shape = 21, color = "black", cex = 4) +
  scale_fill_manual(values = cbPalette[c(7,4)]) +
  scale_color_manual(values = cbPalette[c(7,4)]) +
  xlab("Lamina Layer") +
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
    legend.position = "none",
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.line=element_line()
  )

# ggsave("figures/two.eye.comparison.d13C.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

ggplot(data = bdwf.45.df, aes(x = laminae, y = d15N, color = spp_ID, fill = spp_ID)) +
  geom_line(cex = 1.2, alpha = 0.5) +
  geom_point(shape = 21, color = "black", cex = 4) +
  scale_fill_manual(values = cbPalette[c(7,4)]) +
  scale_color_manual(values = cbPalette[c(7,4)]) +
  xlab("Lamina Layer") +
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
    legend.position = "none",
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.line=element_line()
  )

# ggsave("figures/two.eye.comparison.d15N.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

## Muscle and lamina d15N correlations -------------------------------------

ggplot(correlation.df, aes(x = muscle.d15N, y = lamina.d15N)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, color = "black") +
  geom_point(aes(fill = species, shape = species), size = 5, color = "black") +
  xlab(expression(Muscle~tissue~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  ylab(expression(Outermost~lamina~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  scale_x_continuous(limits = c(8,14), breaks = seq(8,14,1), expand = c(0,0.1)) +
  scale_y_continuous(limits = c(8,14), breaks = seq(8,14,1), expand = c(0,0.1)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2)]) +
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

# ggsave("figures/muscle.vs.lamina.d15N.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

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


## Muscle and lamina d13C correlations -------------------------------------

ggplot(correlation.df, aes(x = muscle.d13C, y = lamina.d13C)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, color = "black") +
  geom_point(aes(fill = species, shape = species), size = 5, color = "black") +
  xlab(expression(Muscle~tissue~italic(delta)^13*C~("\211"~VPDB))) +
  ylab(expression(Outermost~lamina~italic(delta)^13*C~("\211"~VPDB))) +
  scale_x_continuous(limits=c(-25,-21), breaks=seq(-25,-21,1), expand=c(0,0.1)) +
  scale_y_continuous(limits=c(-25,-21), breaks=seq(-25,-21,1), expand=c(0,0.1)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2)]) +
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

# ggsave("figures/muscle.vs.lamina.d13C.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

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


## Muscle and fin d15N correlations ----------------------------------------

ggplot(correlation.df, aes(x = muscle.d15N, y = fin.d15N)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, color = "black") +
  geom_point(aes(fill = species, shape = species), size = 5, color = "black") +
  xlab(expression(Muscle~tissue~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  ylab(expression(Fin~tissue~italic(delta)^15*N~("\211"~atmospheric~N[2]))) +
  scale_x_continuous(limits = c(8,14.2), breaks = seq(8,14,1), expand = c(0,0.1)) +
  scale_y_continuous(limits = c(8,14.2), breaks = seq(8,14,1), expand = c(0,0.1)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2)]) +
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

# ggsave("figures/muscle.vs.fin.d15N.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

muscle.fin.d15N.lm <- lm(fin.d15N ~ muscle.d15N, data = correlation.df)
summary(muscle.fin.d15N.lm) # m = 1.14 (p = <0.01), b = -1.29 (p = 0.37), adj. R^2 = 0.96


## Muscle and fin d13C correlations ----------------------------------------

ggplot(correlation.df, aes(x = muscle.d13C, y = fin.d13C)) +
  geom_abline(slope = 1, intercept = 0, cex = 1, lty = 2, color = "black") +
  geom_point(aes(fill = species, shape = species), size = 5, color = "black") +
  xlab(expression(Muscle~tissue~italic(delta)^13*C~("\211"~VPDB))) +
  ylab(expression(Fin~tissue~italic(delta)^13*C~("\211"~VPDB))) +
  scale_x_continuous(limits=c(-25,-19.5), breaks=seq(-25,-20, 1), expand=c(0,0.1)) +
  scale_y_continuous(limits=c(-25,-19.5), breaks=seq(-25,-20, 1), expand=c(0,0.1)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2)]) +
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

# ggsave("figures/muscle.vs.fin.d13C.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

muscle.fin.d13C.lm <- lm(fin.d13C ~ muscle.d13C, data = correlation.df)
summary(muscle.fin.d13C.lm) # m = 1.28 (p = <0.01), b = 8.16 (p = 0.18), adj. R^2 = 0.91

fin.muscle.adjustment.factor <- mean(correlation.df$fin.d13C - correlation.df$muscle.d13C)
fin.muscle.adjustment.factor # -1.572


## Average d15N per lamina -------------------------------------------------

ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d15N, color = species)) +
  geom_line(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab(expression(Average~italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  xlab("Lamina layer") +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2)]) +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1), expand = c(0,0.2)) +
  scale_y_continuous(limits=c(6,13),breaks=seq(6,13,1), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/average.d15N.vs.lamina.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Average d13C per lamina -------------------------------------------------

ggplot(data = average.delta.laminae, aes(x = laminae, y = mean.d13C)) +
  geom_line(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  ylab(expression(Average~italic(delta)^13*C~("\211"~" VPDB"))) +
  xlab("Lamina layer") +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2)]) +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1), expand = c(0,0.2)) +
  scale_y_continuous(limits=c(-27,-21),breaks=seq(-27,-21,1), expand = c(0,0.2)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size =18),
    axis.text = element_text(size=16, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/average.d13C.vs.lamina.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## Average d15N vs average d13C, by lamina layer ---------------------------

ggplot(data = average.delta.laminae, aes(x = mean.d13C, y = mean.d15N)) +
  geom_path(aes(color = species), cex = 1.2, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 7) +
  geom_text(aes(label=laminae, fontface = "bold"), parse = FALSE) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco"))+
  xlab(expression(italic(delta)^13*C~("\211"~" VPDB"))) +
  ylab(expression(italic(delta)^15*N~("\211"~" atmospheric "~N[2]))) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_color_manual(values=cbPalette[c(1,3,2)]) +
  scale_x_continuous(limits=c(-27,-21),breaks=seq(-27,-21,1), expand = c(0,0.1)) +
  scale_y_continuous(limits=c(6,12.6),breaks=seq(6,12,1), expand = c(0,0.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size =20),
    axis.text = element_text(size=18, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/average.d15N.vs.average.d13C.per.lamina.png", device = "png", dpi = "retina", width = 8, height = 8, units = "in")


## Number of lamina vs fish length -------------------------------------

length.dat <- read.csv("data/sample.LW.data.csv")

length.dat <- length.dat %>%
  unite("spp_ID", species:ID, sep= "_", remove = FALSE) %>%
  inner_join(delta.df, length.dat, by = "spp_ID")

lin.model.dat <- length.dat %>%
  group_by(spp_ID) %>%
  summarize("species" = unique(species.x), "n.laminae" = length(laminae), "length_mm" = unique(length_mm), "weight_kg" = unique(weight_kg))
lin.model.dat
lin.model.dat[3,3] <- 7
lin.model.dat[4,3] <- 9

lm_eqn <- function(df){
  lm <- lm(n.laminae ~ length_mm, data = lin.model.dat);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(lm)[1]), digits = 2),
                        b = format(unname(coef(lm)[2]), digits = 1),
                        r2 = format(summary(lm)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

ggplot(data = lin.model.dat, aes(x = length_mm, y = n.laminae)) +
  geom_smooth(method = "lm", se=FALSE, color="black", fullrange = TRUE) +
  geom_text(x = 375, y = 9.2, label = lm_eqn(df), parse = TRUE, size = 5, angle = 18) +
  geom_point(aes(fill = species, shape = species), color = "black", cex = 5) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_x_continuous(limits=c(215,405),breaks=seq(225,400,25), expand = c(0,0)) +
  scale_y_continuous(limits=c(2.8,11.2),breaks=seq(3,11,1), expand = c(0,0.1)) +
  xlab("Fish length (mm)") +
  ylab("Number of laminae") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size =20),
    axis.text = element_text(size=18, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/length.vs.n.lamina.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")

## C:N ratios --------------------------------------------------------------

cn.df <- bulk.data %>% 
  select(c(species, ID, laminae, Conc.C, Conc.N)) %>% 
  unite("spp_ID", species:ID, sep= "_", remove = FALSE) %>% 
  mutate("C.N.ratio" = Conc.C/Conc.N)

ggplot(data = cn.df, aes(x = spp_ID, y = C.N.ratio)) +
  geom_boxplot() +
  geom_hline(yintercept = 3.5, color = "red", size = 1.2) +
  scale_y_continuous(limits = c(2.85,3.55), breaks = seq(2.9,3.5,0.1)) +
  xlab("Species + ID") +
  ylab(expression(C:N~ratio~("%"[carbon]/"%"[nitrogen]))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.y = element_blank(),
    axis.title = element_text(size=16),
    axis.text.x = element_text(vjust = 0.5, angle = 90),
    axis.text = element_text(size=13, color="black"),
    axis.line=element_line()
  )

ggsave("figures/CN.ratios.lamina.png", device = "png", dpi = "retina", width = 8, height = 4.5, units = "in")



## Two-end-member mixing model ---------------------------------------------

mix.model.df <- delta.df %>%
  group_by(spp_ID, species, ID, laminae) %>%
  summarise(
    "d13C" = d13C,
    "adjusted.d13C" = d13C-2,
    "end1.d13C" = -21.6,
    "end2.d13C" = -27.4,
    "end2.ratio" = ((adjusted.d13C-end1.d13C)/(end2.d13C-end1.d13C)),
    "end1.ratio" = 1-end2.ratio,
    "percent.end1" = 100*end1.ratio,
    "percent.end2" = 100*end2.ratio) %>%
  na.omit()

start.end.mix.model <- mix.model.df %>%
  group_by(spp_ID, species, ID) %>%
  filter(laminae == 0 | laminae == max(laminae)) %>%
  summarise(
    "start.d13C" = d13C[laminae == 0],
    "end.d13C" = d13C[laminae != 0],
    "diff.d13C" = end.d13C - start.d13C,
    "start.percent" = percent.end1[laminae == 0],
    "end.percent" = percent.end1[laminae != 0],
    "diff.percent" = end.percent - start.percent)

mix.model.df %>% 
  group_by(spp_ID, species, ID) %>% 
  filter(laminae == 0) %>% 
ggplot(aes(x = spp_ID, y = d13C, fill = species, shape = species)) +
  geom_hline(yintercept = -27.4, cex = 1.2, color = "black", lty = 2) +
  geom_hline(yintercept = -21.6, cex = 1.2, color = "black", lty = 2) +
  geom_point(cex = 5) +
  xlab("Species ID") +
  ylab(expression(Innermost~lamina~italic(delta)^13*C)) +
  scale_y_continuous(limits = c(-34,-20), breaks = seq(-34,-20,2), expand = c(0,0.2)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text.x = element_text(vjust = 0.5, angle = 90),
    axis.text = element_text(size=10, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/initial.d13C.by.spp_ID.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

mix.model.df %>% 
  group_by(spp_ID, species, ID) %>% 
  filter(laminae == 0) %>% 
ggplot(aes(x = spp_ID, y = percent.end1, fill = species, shape = species)) +
  geom_hline(yintercept = 0, cex = 1.2, color = "black", lty = 2) +
  geom_hline(yintercept = 100, cex = 1.2, color = "black", lty = 2) +
  geom_point(cex = 5) +
  xlab("Species ID") +
  ylab(expression("%"~Carbon[marine])) +
  scale_y_continuous(limits = c(-140,100), breaks = seq(-140,100,20), expand = c(0,10)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text.x = element_text(vjust = 0.5, angle = 90),
    axis.text = element_text(size=10, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/percent.marine.by.spp_ID.1.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

start.end.mix.model %>% 
  select(spp_ID, species, ID, start.percent, end.percent) %>% 
  group_by(spp_ID, species, ID) %>% 
  pivot_longer(cols = c(start.percent, end.percent), names_to = "time", values_to = "percent") %>% 
ggplot(aes(x = spp_ID, y = percent, fill = species, shape = species)) +
  geom_hline(yintercept = 0, cex = 1.2, color = "black", lty = 2) +
  geom_hline(yintercept = 100, cex = 1.2, color = "black", lty = 2) +
  geom_point(aes(alpha = time), cex = 5, show.legend = FALSE) +
  xlab("Species ID") +
  ylab(expression("%"~Carbon[marine])) +
  scale_y_continuous(limits = c(-140,100), breaks = seq(-140,100,20), expand = c(0,10)) +
  scale_alpha_discrete(range = c(1,0.4)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text.x = element_text(vjust = 0.5, angle = 90),
    axis.text = element_text(size=10, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/percent.marine.by.spp_ID.2.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")

start.end.mix.model %>% 
  select(spp_ID, species, ID, start.percent, end.percent) %>% 
  group_by(spp_ID, species, ID) %>% 
  pivot_longer(cols = c(start.percent, end.percent), names_to = "time", values_to = "percent") %>% 
ggplot(aes(x = spp_ID, y = percent, fill = species, shape = species)) +
  geom_hline(yintercept = 0, cex = 1.2, color = "black", lty = 2) +
  geom_hline(yintercept = 100, cex = 1.2, color = "black", lty = 2) +
  geom_point(aes(alpha = time), cex = 5, show.legend = FALSE) +
  geom_line(arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  xlab("Species ID") +
  ylab(expression("%"~Carbon[marine])) +
  scale_y_continuous(limits = c(-140,100), breaks = seq(-140,100,20), expand = c(0,10)) +
  scale_alpha_discrete(range = c(1,0.4)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text.x = element_text(vjust = 0.5, angle = 90),
    axis.text = element_text(size=10, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/percent.marine.by.spp_ID.3.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")


## Delta percent marine vs number of lamina --------------------------------

start.end.mix.model %>% 
  select(spp_ID, species, ID, diff.percent) %>% 
  merge(x = ., y = lin.model.dat) %>% 
ggplot(aes(x = n.laminae, y = diff.percent, fill = species, shape = species)) +
  geom_hline(yintercept = 0, cex = 1.2, color = "black") +
  geom_point(cex = 5) +
  xlab("Total number of laminae") +
  ylab(expression(Delta[outer-inner]~"%"~Carbon[marine])) +
  scale_x_continuous(limits = c(3,11), breaks = seq(3,11,1), expand = c(0,0.3)) +
  scale_y_continuous(limits = c(-20,190), breaks = seq(-20,180,20), expand = c(0,5)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text = element_text(size=10, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/delta.percent.marine.by.n.laminae.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")


## Delta percent marine vs fish length -------------------------------------

start.end.mix.model %>% 
  select(spp_ID, species, ID, diff.percent) %>% 
  merge(x = ., y = lin.model.dat) %>% 
ggplot(aes(x = length_mm, y = diff.percent, fill = species, shape = species)) +
  geom_hline(yintercept = 0, cex = 1.2, color = "black") +
  geom_point(cex = 5) +
  xlab("Fish length (mm)") +
  ylab(expression(Delta[outer-inner]~"%"~Carbon[marine])) +
  scale_x_continuous(limits = c(220,390), breaks = seq(225,375,25), expand = c(0,0)) +
  scale_y_continuous(limits = c(-20,190), breaks = seq(-20,180,20), expand = c(0,5)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text = element_text(size=10, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/delta.percent.marine.vs.length.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")


## Delta percent marine vs fish weight -------------------------------------

start.end.mix.model %>% 
  select(spp_ID, species, ID, diff.percent) %>% 
  merge(x = ., y = lin.model.dat) %>% 
ggplot(aes(x = weight_kg*1000, y = diff.percent, fill = species, shape = species)) +
  geom_hline(yintercept = 0, cex = 1.2, color = "black") +
  geom_point(cex = 5) +
  xlab("Fish weight (g)") +
  ylab(expression(Delta[outer-inner]~"%"~Carbon[marine])) +
  scale_x_continuous(limits = c(100,670), breaks = seq(100,650,50), expand = c(0,0)) +
  scale_y_continuous(limits = c(-20,190), breaks = seq(-20,180,20), expand = c(0,5)) +
  scale_shape_manual(values=c(21,22,24), name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  scale_fill_manual(values=cbPalette[c(1,3,2)], name = "Species", labels = c("Broad whitefish", "Humpback whitefish", "Least cisco")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(color = "gray85"),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text = element_text(size=10, color="black"),
    legend.position = "none",
    axis.line=element_line()
  )

# ggsave("figures/delta.percent.marine.vs.weight.png", device = "png", dpi = "retina", width = 4.5, height = 4.5, units = "in")


# Power analysis ----------------------------------------------------------

## d15N values

between.vars.d15N <- c(rep(NA, 1+length(unique(delta.df$laminae))))
within.vars.d15N <- c(rep(NA, 1+length(unique(delta.df$laminae))))

i = 0
for (i in 0:length(unique(delta.df$laminae))) {
  
  temp.df <- delta.df %>% filter(laminae == i)
  temp.lm <- lm(d15N ~ species, temp.df)
  temp.anova <- anova(temp.lm)
  
  between.vars.d15N[i+1] <- temp.anova$`Mean Sq`[1]
  within.vars.d15N[i+1] <- temp.anova$`Mean Sq`[2]
}

between.vars.d15N
within.vars.d15N

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d15N[1],
                 within.var = within.vars.d15N[1]) # n = 3.466879

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d15N[2],
                 within.var = within.vars.d15N[2]) # n = 2.060118

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d15N[3],
                 within.var = within.vars.d15N[3]) # n = <2

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d15N[4],
                 within.var = within.vars.d15N[4]) # n = <2

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d15N[5],
                 within.var = within.vars.d15N[5]) # n = <2

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d15N[6],
                 within.var = within.vars.d15N[6]) # n = <2

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d15N[7],
                 within.var = within.vars.d15N[7]) # n = <2

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8, 
                 between.var = between.vars.d15N[8],
                 within.var = within.vars.d15N[8]) # n = <2

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8, 
                 between.var = between.vars.d15N[9],
                 within.var = within.vars.d15N[9]) # n = 4.322733


## d13C values

between.vars.d13C <- c(rep(NA, 1+length(unique(delta.df$laminae))))
within.vars.d13C <- c(rep(NA, 1+length(unique(delta.df$laminae))))

i = 0
for (i in 0:length(unique(delta.df$laminae))) {
  
  temp.df <- delta.df %>% filter(laminae == i)
  temp.lm <- lm(d13C ~ species, temp.df)
  temp.anova <- anova(temp.lm)
  
  between.vars.d13C[i+1] <- temp.anova$`Mean Sq`[1]
  within.vars.d13C[i+1] <- temp.anova$`Mean Sq`[2]
}

between.vars.d13C
within.vars.d13C

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d13C[1],
                 within.var = within.vars.d13C[1]) # n = 2.558315

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d13C[2],
                 within.var = within.vars.d13C[2]) # n = 2.903184

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d13C[3],
                 within.var = within.vars.d13C[3]) # n = 2.678831

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d13C[4],
                 within.var = within.vars.d13C[4]) # 2.921139

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d13C[5],
                 within.var = within.vars.d13C[5]) # n = <2

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d13C[6],
                 within.var = within.vars.d13C[6]) # n = 2.006852

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8,
                 between.var = between.vars.d13C[7],
                 within.var = within.vars.d13C[7]) # n = 4.371279

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8, 
                 between.var = between.vars.d13C[8],
                 within.var = within.vars.d13C[8]) # n = 16.31604

power.anova.test(groups = 4, sig.level = 0.05, power = 0.8, 
                 between.var = between.vars.d13C[9],
                 within.var = within.vars.d13C[9]) # n = 14.18908

## Back-of-the-envelope sample-cost calculations:

delta.df %>% 
  group_by(species) %>% 
  summarise("count" = length(unique(ID)))
adult.fish.summer1 <- 8*1*10 + 3*3*10 #170
adult.fish.summer2 <- 8*4*10 #320
juvenile.fish.summer2 <- 5*4*5 #100

total.lam.samples.bulk <- sum(adult.fish.summer1, adult.fish.summer2, juvenile.fish.summer2) # 590
total.bulk.costs.left <- 14*total.lam.samples.bulk # $8,260
