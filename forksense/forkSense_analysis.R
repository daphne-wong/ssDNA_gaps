#### Load libraries ####
library(ggplot2)

#### Set Directories and Sample ####

setwd("/project2/tp_612_1653/for_Tirzah/dnascent_petrina/forksense/minLength_100")
getwd()

out_dir <- paste0(getwd(),"/plots/")
out_dir

#### Read files ####
in_dir <- paste0(getwd(),"/parsed_output/")
in_dir

files <- list.files(path=in_dir, pattern = "^ALL_samples_.*\\.tsv$")
files
var_names <- sub("^ALL_samples_", "", files)
var_names <- sub("\\.tsv$", "", var_names)
var_names

# Read each file and assign it to a variable of the corresponding name
for (i in seq_along(files)) {
  df <- read.delim(paste0(in_dir,files[i]), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  assign(var_names[i], df)
}



#### Check segment size ####

EdU_segments$segment_size <- EdU_segments$EdU_end - EdU_segments$EdU_start
summary(EdU_segments$segment_size)
EdU_1000 <- sum(EdU_segments$segment_size < 1000)
EdU_2000 <- sum(EdU_segments$segment_size > 2000)

message("Number of EdU Segments < 1000 bp:   ", EdU_1000," out of ", nrow(EdU_segments))
message("Number of EdU Segments > 2000 bp:   ", EdU_2000, " out of ", nrow(EdU_segments))

BrdU_segments$segment_size <- BrdU_segments$BrdU_end - BrdU_segments$BrdU_start
summary(BrdU_segments$segment_size)
BrdU_1000 <- sum(BrdU_segments < 1000)
BrdU_2000 <- sum(BrdU_segments > 2000)

message("Number of BrdU Segments < 1000 bp:   ", BrdU_1000," out of ", nrow(BrdU_segments))
message("Number of BrdU Segments > 2000 bp:   ", BrdU_2000, " out of ", nrow(BrdU_segments))



#### Finding flanking segments ####



# Find how 










#### Plotting ####


