################################################################################
# Section 4.3 Case study
# Prepare moss data
################################################################################
rm(list = ls())

library(readr)
library(ggplot2)
library(sf)
library(dplyr)

inpath  <- "../data/"
outpath <- "../outputs/real-bh-moss/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# Load survey grid locations
#-------------------------------------------------------------------------------
survey_locs <- read_csv(paste0(outpath, "survey_locs.csv"))
survey_locs$Eastings <- survey_locs$Eastings*100 + 500000
survey_locs$Northings <- survey_locs$Northings*100 + 2600000

survey_eastnorth_sf <- st_as_sf(survey_locs, coords = c("Eastings", "Northings"),
                                crs = "+proj=utm +zone=47 +south +datum=WGS84 +units=m +no_defs")
survey_lonlat_sf <- st_transform(survey_eastnorth_sf, 
                                 crs = "+proj=longlat +datum=WGS84 +no_defs")
survey_lonlat <- st_coordinates(survey_lonlat_sf)
survey_locs$lon <- survey_lonlat[, 1]
survey_locs$lat <- survey_lonlat[, 2]


#-------------------------------------------------------------------------------
# Load moss distribution data
#-------------------------------------------------------------------------------
bh_moss <- read_csv(paste0(inpath, "moss/bh_moss.csv"))
bh_moss <- bh_moss %>% arrange(Eastings, Northings) %>% filter(!is.na(Northings))
bh_moss <- bh_moss[!duplicated(bh_moss[, 1:2]), ]

bh_moss$moss <- ifelse(is.na(bh_moss$Moss), 0, 1)
bh_moss <- bh_moss[, c(1,2,4)]


#-------------------------------------------------------------------------------
# Combine data
#-------------------------------------------------------------------------------
bh_moss$Eastings <- bh_moss$Eastings*100 + 500000
bh_moss$Northings <- bh_moss$Northings*100 + 2600000
survey_moss <- left_join(survey_locs, bh_moss, by = c("Eastings", "Northings"))
survey_moss$moss[is.na(survey_moss$moss)] <- 0

save(survey_moss, file = paste0(outpath, "survey_moss.rdata"))


#-------------------------------------------------------------------------------
# Load Bunger Hills map
#-------------------------------------------------------------------------------
bungerhills <- st_read(paste0(inpath, "bh/BH_mask.shp"))
map_trunc <- c(564500, 596500, 2639500, 2660500)
bh_map <- bungerhills %>% st_crop(c(xmin = 564500, ymin = 2639500, xmax = 596500, ymax = 2660500))
  

#-------------------------------------------------------------------------------
# Plot moss distribution
#-------------------------------------------------------------------------------
survey_moss$moss2 <- ifelse(survey_moss$moss == 1, "presence", "absence")
png(paste0(outpath, "bh_survey_moss", ".png"), width=1000, height=600, pointsize=20, bg = "transparent")
ggplot(data = bh_map) + geom_sf() + coord_sf(datum = st_crs(bh_map)) + 
  geom_point(aes(x = Eastings, y = Northings, col = moss2, shape = moss2), data = survey_moss, size = 2) +
  # geom_point(aes(x = Eastings, y = Northings, shape = moss2), data = survey_moss, size = 3) +
  scale_shape_manual(values = c(1, 2)) + 
  xlab("Eastings") + ylab("Northings") + 
  labs(col = "moss (all speices)", shape = "moss (all speices)") +
  # labs(shape = "moss") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm"), angle = 90),
        legend.key.height= unit(1.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
dev.off()


png(paste0(outpath, "survey_grid_cells", ".png"), width=1000, height=600, pointsize=20, bg = "transparent")
ggplot(data = bh_map) + geom_sf() + coord_sf(datum = st_crs(bh_map)) + 
  geom_point(aes(x = Eastings, y = Northings), data = survey_moss, size = 3) +
  scale_x_continuous(breaks = c(570000, 580000, 590000)) +
  scale_y_continuous(breaks = c(2640000, 2650000, 2660000)) +
  scale_shape_manual(values = c(1, 19)) + 
  xlab("Eastings") + ylab("Northings") + 
  theme(axis.text = element_text(size = 28, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 40, b = 0, l = 0)),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(margin = unit(c(0.5, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.5, 0, 0.2), "cm"), angle = 90, hjust = 0.5, vjust = 0.5),
        plot.margin = margin(t = 55, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.height= unit(1.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
dev.off()


