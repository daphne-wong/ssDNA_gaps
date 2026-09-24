#### Setup ####
getwd()
setwd("/Users/daphnewong/Documents/USC/ssDNA_gaps")
getwd()
library(ggplot2)
library(dplyr)
library(scales)


gaps_dir <- paste0(getwd(),"/gaps_python/20260921_gap_info")
gaps_dir
sample_list <- c("UWB_Cisp0", "UWB_Cisp1", "OV3_Cisp0", "OV3_Cisp1_masked")

#### Read data into dfs ####

all.files <- list.files(path=gaps_dir, pattern=".*\\.tsv$", full.names=T)
all.files

split_sample_name <- function(sample_names) {
  parts <- strsplit(sample_names, "_")
  cell_line <- sapply(parts, `[`, 1)
  treatment <- sapply(parts, function(x) paste(x[-1], collapse = "_"))
  treatment <- sub("_masked$", "", treatment)
  data.frame(cell_line = cell_line, treatment = treatment, sample = sample_names,
             stringsAsFactors = FALSE)
}

data <- list()

for (f in all.files) {
  base_name <- sub("\\.tsv$", "", basename(f))
  type <- sub("^.*\\.", "", base_name)
  sample <- sub("\\.[^.]+$", "", base_name)
  
  if (is.null(data[[type]])) data[[type]] <- list()
  data[[type]][[sample]] <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

for (type in names(data)) {
  dfs <- data[[type]]
  meta <- split_sample_name(names(dfs))
  
  assign(type, do.call(rbind, lapply(seq_along(dfs), function(i) {
    cbind(meta[i, , drop = FALSE], dfs[[i]], row.names = NULL)
  })))
}

rm(dfs, data)




#### Analysis ####


# Number of gaps per sample
per_sample_gaps <- 


gaps$chrom <- factor(gaps$chrom, levels = paste0("chr", c(1:22, "X", "Y")))

per_chrom$chrom <- factor(
  per_chrom$chrom,
  levels = paste0("chr", c(1:22, "X", "Y"))
)

# GRCh38 chromosome lengths
chr_lengths <- data.frame(
  chrom = paste0("chr", c(1:22, "X", "Y")),
  chr_length_bp = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468, 156040895, 57227415
  )
)

per_chrom$chr_length_bp <- chr_lengths$chr_length_bp[match(as.character(per_chrom$chrom), chr_lengths$chrom)]
per_chrom$gaps_per_chr_Mb <- per_chrom$n_gaps / (per_chrom$chr_length_bp / 1e6)
subset(per_chrom, is.na(chr_length_bp))

chrom_rates <- per_chrom %>%
  group_by(sample, chrom) %>%
  summarise(gap_rate = sum(n_gaps) / sum(interrogated_Mb),
            .groups = "drop") %>%
  mutate(chrom = factor(chrom, levels = paste0("chr", c(1:22, "X", "Y"))))

per_read2 <- per_read %>%
  mutate(read_span = ref_end - ref_start + 1,
         has_gap = n_gaps > 0,
         gaps_per_Mb = n_gaps / informative_bp * 1e6)



# Per chromosome summary
summary_per_chr <- per_chrom %>%
  group_by(cell_line, treatment, sample, chrom) %>%
  summarise(across(where(is.numeric),~ mean(.x, na.rm = TRUE)),
            .groups = "drop")
summary_per_chr



# Percentage of reads that contain a gap
gap_fraction <- per_read2 %>%
  group_by(cell_line, treatment, sample) %>%
  summarise(prop_with_gap = mean(has_gap),
            n_reads = n(),
            .groups = "drop")
gap_fraction


# Gap rate by sample
sample_rates <- per_read2 %>%
  group_by(cell_line, treatment, sample) %>%
  summarise(gap_rate = sum(n_gaps) / sum(informative_bp) * 1e6,
            .groups = "drop")


### Plots ####

# Distribution of gap lengths for each sample
ggplot(gaps, aes(x = length, fill = treatment, color = treatment)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~ sample, ncol = 2, scales = "free_y") +
  labs(x = "Gap length (bp)", y = "Density", title="Distribution of gap lengths", fill = "Treatment", color = "Treatment") +
  theme_classic()

ggplot(subset(gaps, length > 0), aes(x = length, fill = treatment, color = treatment)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_x_log10(labels = label_number()) +
  facet_wrap(~ sample, ncol = 2, scales = "free_y") +
  labs(x = "Gap length (bp, log scale)", y = "Density", title="Distribution of gap lengths (log-scaled)", fill = "Treatment", color = "Treatment") +
  theme_classic()

# N gaps by sample
ggplot(gaps, aes(x = length, fill = treatment)) +
  geom_histogram(binwidth = 50, color = "white", boundary = 0) +
  facet_wrap(~ sample, ncol = 2, scales = "free_y") +
  labs(x = "Gap length (bp)", y = "Number of gaps", fill = "Treatment") +
  theme_classic()


# N_gaps by chromosome
ggplot(per_chrom, aes(x = chrom, y = n_gaps)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ sample, ncol = 2) +
  labs(x = "Chromosome", y = "Number of gaps", title="n_gaps by Chromosome") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# N_gaps by chromosome, normalized by physical chromosome length
ggplot(per_chrom, aes(x = chrom, y = gaps_per_chr_Mb)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ sample, ncol = 2) +
  labs(x = "Chromosome", y = "Observed gaps per chromosome Mb", title = "Gap density normalized by chromosome length") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# gaps_per_Mb by chromosome (normalised by interrogated_Mb length)
ggplot(per_chrom, aes(x = chrom, y = gaps_per_Mb)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ sample, ncol = 2) +
  labs(x = "Chromosome", y = "Gaps_per_Mb", title="Gaps_per_Mb by Chromosome") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Chromosome-by-sample gap-rate heat map
# chrY has highest density of gaps per Mb (although it is also the smallest)

ggplot(chrom_rates, aes(x = chrom, y = sample, fill = gap_rate)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Gaps\nper Mb") +
  labs(x = "Chromosome", y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Cumulative proportion of gaps per Mb per chromosome
# For example, >90% of chromosomes in Cisp1 has 
ggplot(per_chrom, aes(x = gaps_per_Mb, color = treatment)) +
  stat_ecdf(linewidth = 1) +
  facet_wrap(~ cell_line) +
  labs(x = "Gaps per Mb", y = "Cumulative proportion of chromosomes", title="Cumulative proportion of gaps per Mb per chromosome", color = "Treatment") +
  theme_classic() +
  theme(legend.position = "top")

# Distribution of calls_per_kb across gaps (normalized by individual gap length)
# Within each sample and chromosome, how densely packed are the calls supporting individual gap regions?
# How densely packed are the calls supporting each gap region?
ggplot(gaps, aes(x = chrom, y = calls_per_kb, fill = treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.alpha = 0.15, width = 0.7) +
  facet_wrap(~ cell_line, ncol = 1) +
  labs(x = "Chromosome", y = "Call density (calls/kb of gap)", title = "Call density within gaps by chromosome", fill = "Treatment") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(gaps, aes(x = chrom, y = calls_per_kb, fill = treatment)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.6) + 
  facet_wrap(~ sample) +
  labs(x = "Chromosome", y = "Calls per kb", fill = "Treatment") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







# Distribution of gap rates among reads
ggplot(filter(per_read2, has_gap), aes(x = gaps_per_Mb, color = treatment)) +
  stat_ecdf(linewidth = 1) +
  scale_x_log10() +
  facet_wrap(~ cell_line) +
  labs(x = "Gaps per Mb among gap-containing reads", y = "Cumulative proportion", color = "Treatment") +
  theme_classic()




