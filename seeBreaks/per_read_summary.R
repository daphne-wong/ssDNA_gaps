options(scipen = 999)

per_read_files <- list.files(
  path = ".",
  pattern = "_forkSense_per_read\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


#### Function to calculate statistics ####

calculate_stats <- function(x) {
  x <- x[!is.na(x)]
  
  quartiles <- quantile(
    x,
    probs = c(0, 0.25, 0.50, 0.75, 1),
    names = FALSE
  )
  
  data.frame(
    n_reads = length(x),
    n_positive_reads = sum(x > 0),
    total_labelled_bp = sum(x),
    mean_bp = mean(x),
    sd_bp = if (length(x) > 1) sd(x) else NA_real_,
    minimum_bp = quartiles[1],
    q1_bp = quartiles[2],
    q2_median_bp = quartiles[3],
    q3_bp = quartiles[4],
    q4_maximum_bp = quartiles[5]
  )
}


#### Sample-level statistics ####

sample_results <- lapply(per_read_files, function(file) {
  
  dat <- read.delim(file, check.names = FALSE)
  
  if (nrow(dat) == 0) {
    return(NULL)
  }
  
  sample_name <- unique(dat$sample)[1]
  
  rbind(
    cbind(
      sample = sample_name,
      analogue = "EdU",
      calculate_stats(dat$EdU_bp)
    ),
    cbind(
      sample = sample_name,
      analogue = "BrdU",
      calculate_stats(dat$BrdU_bp)
    ),
    cbind(
      sample = sample_name,
      analogue = "EdU+BrdU",
      calculate_stats(dat$total_analogue_bp)
    )
  )
})

summary_table <- do.call(rbind, sample_results)
rownames(summary_table) <- NULL


#### Per-chromosome statistics ####

chromosome_results <- lapply(per_read_files, function(file) {
  
  dat <- read.delim(file, check.names = FALSE)
  
  if (nrow(dat) == 0) {
    return(NULL)
  }
  
  sample_name <- unique(dat$sample)[1]
  chromosome_groups <- split(dat, dat$chromosome)
  
  do.call(
    rbind,
    lapply(names(chromosome_groups), function(chr) {
      
      chr_dat <- chromosome_groups[[chr]]
      
      rbind(
        cbind(
          sample = sample_name,
          chromosome = chr,
          analogue = "EdU",
          calculate_stats(chr_dat$EdU_bp)
        ),
        cbind(
          sample = sample_name,
          chromosome = chr,
          analogue = "BrdU",
          calculate_stats(chr_dat$BrdU_bp)
        ),
        cbind(
          sample = sample_name,
          chromosome = chr,
          analogue = "EdU+BrdU",
          calculate_stats(chr_dat$total_analogue_bp)
        )
      )
    })
  )
})

chromosome_summary_table <- do.call(rbind, chromosome_results)
rownames(chromosome_summary_table) <- NULL


#### Format numeric results ####

format_summary <- function(dat) {
  
  dat$total_labelled_bp <- round(dat$total_labelled_bp)
  
  decimal_columns <- c(
    "mean_bp",
    "sd_bp",
    "minimum_bp",
    "q1_bp",
    "q2_median_bp",
    "q3_bp",
    "q4_maximum_bp"
  )
  
  dat[decimal_columns] <- lapply(
    dat[decimal_columns],
    round,
    digits = 2
  )
  
  dat
}

summary_table <- format_summary(summary_table)
chromosome_summary_table <- format_summary(chromosome_summary_table)


#### Save TSV files ####

write.table(
  summary_table,
  file = "per_sample_analogue_bp_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  chromosome_summary_table,
  file = "per_chr_analogue_bp_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)


#### Display results ####

summary_table
chromosome_summary_table