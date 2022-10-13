# =====================================================================
# Nitrogen amino acid compound-specific SIA from eye lenses
# Jonah Bacon
# 12 October 2022
# =====================================================================


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggsci)
library(ggrepel)

# Load data ---------------------------------------------------------------

csia.dat <- read.csv("data/CSIA_data_1.csv", header = TRUE, skip = 1)

temp.dat1 <- filter(csia.dat, CSIA_ID != is.na(CSIA_ID))

temp.dat2 <- temp.dat1 %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            Ala = mean(Ala, na.rm = T), Gly = mean(Gly, na.rm = T), Val = mean(Val, na.rm = T), Leu = mean(Leu, na.rm = T), iLeu = mean(iLeu, na.rm = T),
            nLeu = mean(nLeu, na.rm = T), Thr = mean(Thr, na.rm = T), Ser = mean(Ser, na.rm = T), Asp = mean(Asp, na.rm = T), Pro = mean(Pro, na.rm = T),
            Glu = mean(Glu, na.rm = T), Met = mean(Met, na.rm = T), Phe = mean(Phe, na.rm = T), Caf = mean(Caf, na.rm = T), Lys = mean(Lys, na.rm = T))

weighted.means <- temp.dat1 %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            del.15N.Gly = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T)) / sum(1/sd(Gly, na.rm = T)), 
            del.15N.Phe = sum(mean(Phe, na.rm = T)/sd(Phe, na.rm = T)) / sum(1/sd(Phe, na.rm = T)), 
            del.15N.Lys = sum(mean(Lys, na.rm = T)/sd(Lys, na.rm = T)) / sum(1/sd(Lys, na.rm = T)),
            del.15N.Ala = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T)) / sum(1/sd(Ala, na.rm = T)), 
            del.15N.Leu = sum(mean(Leu, na.rm = T)/sd(Leu, na.rm = T)) / sum(1/sd(Leu, na.rm = T)), 
            del.15N.Glu = sum(mean(Glu, na.rm = T)/sd(Glu, na.rm = T)) / sum(1/sd(Glu, na.rm = T)))

beta = 3.6
TEF = 5.7

ARCS.weighted.means <- temp.dat1 %>% 
  filter(Species == "ARCS") %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.nom = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.nom = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denom = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denom = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.nom/source.denom,
            trophic = trophic.nom/trophic.denom,
            TP = ((trophic - source - beta)/TEF) + 1)

BDWF.weighted.means <- temp.dat1 %>% 
  filter(Species == "BDWF") %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.nom = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.nom = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denom = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denom = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.nom/source.denom,
            trophic = trophic.nom/trophic.denom,
            TP = ((trophic - source - beta)/TEF) + 1)

LSCS.weighted.means <- temp.dat1 %>% 
  filter(Species == "LSCS") %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.nom = sum(mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.nom = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denom = sum(1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denom = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.nom/source.denom,
            trophic = trophic.nom/trophic.denom,
            TP = ((trophic - source - beta)/TEF) + 1)

HBWF1.weighted.means <- temp.dat1 %>% 
  filter(Species == "HBWF" & CSIA_ID == 8) %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.nom = sum(mean(Gly, na.rm = T)/sd(Gly, na.rm = T), mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.nom = sum(mean(Ala, na.rm = T)/sd(Ala, na.rm = T), mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denom = sum(1/sd(Gly, na.rm = T), 1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denom = sum(1/sd(Ala, na.rm = T), 1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.nom/source.denom,
            trophic = trophic.nom/trophic.denom,
            TP = ((trophic - source - beta)/TEF) + 1)

HBWF2.weighted.means <- temp.dat1 %>% 
  filter(Species == "HBWF" & CSIA_ID == 10) %>% 
  group_by(CSIA_ID, Sample_ID) %>% 
  summarise(Species = unique(Species), ID = unique(ID), spp_ID = unique(spp_ID), Layer = unique(Layer), Sample_wt_mg = unique(Sample_wt_mg),
            source.nom = sum(mean(Phe, na.rm = T)/sd(Phe, na.rm = T), mean(Lys, na.rm = T)/sd(Lys, na.rm = T)),
            trophic.nom = sum(mean(Leu, na.rm = T)/sd(Leu, na.rm = T), mean(Glu, na.rm = T)/sd(Glu, na.rm = T)),
            source.denom = sum(1/sd(Phe, na.rm = T), 1/sd(Lys, na.rm = T)),
            trophic.denom = sum(1/sd(Leu, na.rm = T), 1/sd(Glu, na.rm = T)),
            source = source.nom/source.denom,
            trophic = trophic.nom/trophic.denom,
            TP = ((trophic - source - beta)/TEF) + 1)


combined.csia.weighted.means.df <- rbind(
  ARCS.weighted.means,
  LSCS.weighted.means,
  BDWF.weighted.means,
  HBWF1.weighted.means,
  HBWF2.weighted.means
)
combined.csia.weighted.means.df$tempvar <- c(2,1,2,1,3,1,2,2,1)

ggplot(combined.csia.weighted.means.df, aes(x = tempvar, y = TP, color = Species, fill = Species, label = round(TP,2))) +
  scale_color_jco() +
  scale_fill_jco() +
  geom_line(lwd = 2) +
  geom_point(pch = 21, size = 4, color = "black") +
  geom_label_repel(fill = "white", color = "black", box.padding = 1) +
  xlab("Relative Time Point") +
  ylab("Trophic Position") +
  scale_x_continuous(limits = c(0.8,3.2), breaks = seq(1,3,1), expand = c(0,0)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "gray70"),
    text = element_text(size = 14, face = "bold")
  )
