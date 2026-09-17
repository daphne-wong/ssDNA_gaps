#### Set Directories and Sample ####

library(ggplot2)
library(dplyr)

samples <- c("OV3_Cisp0","OV3_Cisp1_masked","UWB_Cisp0","UWB_Cisp1")

setwd("/project2/tp_612_1653/for_Tirzah/dnascent_petrina/seeBreaks/")
getwd()

break_dir <- "/project2/tp_612_1653/for_Tirzah/dnascent_petrina/seeBreaks/"
break_dir


#### Parse seeBreaks output ####

# tail -n +30017 UWB_Cisp1_seeBreaks.dnascent > UWB_Cisp1_ObservedReadEndFractions.dnascent
# awk 'NR >= 16 && NR <= 30015' OV3_Cisp1_masked_seeBreaks.dnascent > OV3_Cisp1_masked_ExpectedReadEndFractions.dnascent


# Read Observed Read End Fractions
samples_obsREF <- setNames(
  lapply(samples, function(sample) {read.table(paste0(break_dir, sample, "_ObservedReadEndFractions.dnascent"),
                                               col.names = "observed_read_end_fraction")}),
  samples
)

# Read Expected Read End Fractions
samples_expREF <- setNames(
  lapply(samples, function(sample) {read.table(paste0(break_dir, sample, "_ExpectedReadEndFractions.dnascent"),
                                               col.names = "expected_read_end_fraction")}),
  samples
)


lapply(samples_obsREF, function(df) {summary(df)})
lapply(samples_expREF, function(df) {summary(df)})




#### Plot read-end distributions for each df ####

# Combine all obs data into 1 df
obs_plot_df <- bind_rows(samples_obsREF, .id = "sample") %>%
  transmute(sample, 
            distribution = "Observed",
            read_end_fraction = observed_read_end_fraction)


# Combine all exp data into 1 df
exp_plot_df <- bind_rows(samples_expREF, .id = "sample") %>%
  transmute(sample,
            distribution = "Expected",
            read_end_fraction = expected_read_end_fraction)


# Combine both obs and exp into 1 df
plot_df <- bind_rows(obs_plot_df, exp_plot_df)


# Plot with ggplot2

#### Density plot ####
ggplot(plot_df, aes(x = read_end_fraction,
                    color = distribution,
                    fill = distribution)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  facet_wrap(~sample, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Read-end fraction",
       y = "Density",
       color = NULL,
       fill = NULL) +
  theme_classic() +
  theme(legend.position = "top")


#### Boxplot ####
ggplot(plot_df,aes(x = distribution, 
                   y = read_end_fraction,
                   fill = distribution)) +
  geom_boxplot(width = 0.65,outlier.alpha = 0.25) +
  facet_wrap(~sample) +
  coord_cartesian(ylim = c(0, 0.75)) +
  labs(x = NULL,
       y = "Read-end fraction",
       fill = NULL) +
  theme_classic() +
  theme(legend.position = "none")
  
