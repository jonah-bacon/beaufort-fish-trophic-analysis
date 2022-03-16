#'---
#' title: "Trophic niche analysis of four Arctic whitefishes"
#' author: "Jonah Bacon"
#' output: 
#'    html_document:
#'       toc: true
#'       theme: default
#'       toc_depth: 2
#'       toc_float:
#'           collapsed: false
#'---

# header ------------------------------------------------------------------

#' # Preface
#' 
#' This code will analyze gut content data and stable isotope data of four Corregonidae species. Fish were sampled in the summer of 2021 in the nearshore environment of the Beaufort Sea near Prudhoe Bay, AK. This data aims to determine the degree of trophic overlap and niche partitioning currently occurring within these four species. 

# load libraries ----------------------------------------------------------
#' ### Load libraries:
#' 
  library(ggplot2)
  library(data.table)


# load data ---------------------------------------------------------------
#' # Load and Review Data
#' ## Load data files

  # load length:weight data set
  dat <- fread(file = "data/data.csv")

# clean up & explore data -------------------------------------------------

#' ## Clean up L:W data
  names(dat)
  str(dat)

  # make species and net_ID column values factors
  dat$species <- as.factor(dat$species)  
  dat$net_ID <- as.factor(dat$net_ID)  

  # change date formats
  dat[ , date_collected := as.IDate(date_collected, format = "%d-%b-%Y"), ]
  dat[ , date_processed := as.IDate(date_processed, format = "%d-%b-%Y"), ]
  
  # make length_mm column values numeric
  dat$length_mm <- as.numeric(dat$length_mm)  
  
#' ## Visualize L:W relationships

  ggplot(dat, aes(x = length_mm, y = weight_kg)) +
    geom_point() +
    geom_smooth(method = lm) +
    facet_grid(vars(group = species)) +
    xlab("Length (mm)") +
    ylab("Weight (kg)")
  
  boxplot(dat$length_mm ~ dat$species, xlab = "Species", ylab = "Length (mm)")
  
  ggplot(dat) +
    geom_density(aes(length_mm, col = species)) +
    xlab("Length (mm)") +
    ylab("Density")

  boxplot(dat$weight_kg ~ dat$species, xlab = "Species", ylab = "Weight (kg)")
  
  ggplot(dat) +
    geom_density(aes(weight_kg, col = species)) +
    xlab("Weight (kg)") +
    ylab("Density")    
  

  