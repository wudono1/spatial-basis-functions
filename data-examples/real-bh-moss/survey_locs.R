################################################################################
# Section 4.3 Case study
# Prepare moss data
################################################################################
rm(list = ls())

library(dplyr)
library(readr)

inpath  <- "../data/"
outpath <- "../outputs/real-bh-moss/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# Leishman et al. (2020) Figure 2 survey locations
#-------------------------------------------------------------------------------
col01 <- cbind(665, c(455, 465, 475, 485, 505))
col02 <- cbind(675, seq(445, 505, by = 10))
col03 <- cbind(685, seq(435, 505, by = 10))
col04 <- cbind(695, seq(435, 505, by = 10))
col05 <- cbind(705, c(seq(435, 505, by = 10), 525, 535))
col06 <- cbind(715, c(seq(435, 515, by = 10), 535))
col07 <- cbind(725, seq(435, 525, by = 10))
col08 <- cbind(735, seq(415, 545, by = 10))
col09 <- cbind(745, c(415, seq(435, 565, by = 10)))
col10 <- cbind(755, seq(425, 575, by = 10))
col11 <- cbind(765, seq(425, 565, by = 10))
col12 <- cbind(775, seq(425, 575, by = 10))
col13 <- cbind(785, c(seq(415, 515, by = 10), seq(535, 575, by = 10)))
col14 <- cbind(795, c(seq(415, 495, by = 10), 545))
col15 <- cbind(805, c(415, 425, 435, seq(455, 515, by = 10)))
col16 <- cbind(815, c(415, 425, 435, seq(455, 525, by = 10)))
col17 <- cbind(825, seq(415, 535, by = 10))
col18 <- cbind(835, seq(415, 545, by = 10))
col19 <- cbind(845, seq(415, 515, by = 10))
col20 <- cbind(855, seq(415, 525, by = 10))
col21 <- cbind(865, seq(425, 515, by = 10))
col22 <- cbind(875, c(455, 465, 475, 485, 495, 515, 575, 585))
col23 <- cbind(885, c(455, 465, 475, 485, 495, 505, 515, 585))
col24 <- cbind(895, c(465, 475, 485, 495, 505, 575))
col25 <- cbind(905, c(485, 505))
col26 <- cbind(915, 505)
col27 <- cbind(925, c(495, 515, 525))

fig_locs <- rbind(col01, col02, col03, col04, col05,
                  col06, col07, col08, col09, col10,
                  col11, col12, col13, col14, col15,
                  col16, col17, col18, col19, col20,
                  col21, col22, col23, col24, col25,
                  col26, col27)
colnames(fig_locs) <- c("Eastings", "Northings")


#-------------------------------------------------------------------------------
# Bunger Hills data set locations
#-------------------------------------------------------------------------------
bh_moss <- read_csv(paste0(inpath, "moss/bh_moss.csv"))
bh_moss <- bh_moss %>% arrange(Eastings, Northings) %>% filter(!is.na(Northings)) 
bh_moss <- bh_moss %>% distinct(Eastings, Northings)
bh_moss_locs <- bh_moss[, 1:2]


#-------------------------------------------------------------------------------
# Combine
#-------------------------------------------------------------------------------
survey_locs <- dplyr::union(as.data.frame(fig_locs), bh_moss_locs)
write_csv(survey_locs, file = paste0(outpath, "survey_locs.csv"))
