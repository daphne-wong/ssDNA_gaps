#### Set Directories and Sample ####

library(ggplot2)
library(dplyr)

samples <- c("OV3_Cisp0","OV3_Cisp1_masked","UWB_Cisp0","UWB_Cisp1")

setwd("/project2/tp_612_1653/for_Tirzah/dnascent_petrina/meIODy/")
getwd()

melody_dir <- "/project2/tp_612_1653/for_Tirzah/dnascent_petrina/meIODy/"
melody_dir


#### Parse seeBreaks output ####

# tail -n +30017 UWB_Cisp1_seeBreaks.dnascent > UWB_Cisp1_ObservedReadEndFractions.dnascent
# awk 'NR >= 16 && NR <= 30015' OV3_Cisp1_masked_seeBreaks.dnascent > OV3_Cisp1_masked_ExpectedReadEndFractions.dnascent


# Read Observed Read End Fractions
samples_obsREF <- setNames(
  lapply(samples, function(sample) {read.table(paste0(melody_dir, sample, "_ObservedReadEndFractions.dnascent"),
                                               col.names = "observed_read_end_fraction")}),
  samples
)

# Read Expected Read End Fractions
samples_expREF <- setNames(
  lapply(samples, function(sample) {read.table(paste0(melody_dir, sample, "_ExpectedReadEndFractions.dnascent"),
                                               col.names = "expected_read_end_fraction")}),
  samples
)


lapply(samples_obsREF, function(df) {summary(df)})
lapply(samples_expREF, function(df) {summary(df)})
