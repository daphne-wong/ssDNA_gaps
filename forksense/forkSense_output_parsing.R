#### Set Directories and Sample ####

samples <- c("OV3_Cisp0","OV3_Cisp1_masked","UWB_Cisp0","UWB_Cisp1")

setwd("/project2/tp_612_1653/for_Tirzah/dnascent_petrina/forksense/minLength_100")
getwd()

in_dir <- paste0(getwd(),"/")

out_dir <- paste0(in_dir,"/parsed_output/")
out_dir


# Setup for collating
df_names <- c("BrdU_segments", "EdU_segments", "leftForks", "rightForks", "allForks", "origins", "terminations", "leftStress", "rightStress", "allStress")
collated <- setNames(vector("list", length(df_names)), df_names)

#### Parsing Loop ####

for (sample in samples) {
  sample_dir <- paste0(in_dir,sample,"/")
  sample_dir
  
  #### BrdU and EdU ####
  BrdU_segments <- read.table(paste0(sample_dir,"BrdU_DNAscent_forkSense.bed"))
  EdU_segments <- read.table(paste0(sample_dir,"EdU_DNAscent_forkSense.bed"))
  
  colnames(BrdU_segments) <- c("chromosome","BrdU_start","BrdU_end","readID","read_start","read_end","strand")
  colnames(EdU_segments) <- c("chromosome","EdU_start","EdU_end","readID","read_start","read_end","strand")
  
  
  #### Forks ####
  leftForks <- read.table(paste0(sample_dir,"leftForks_DNAscent_forkSense.bed"))
  rightForks <- read.table(paste0(sample_dir,"rightForks_DNAscent_forkSense.bed"))
  
  colnames(leftForks) <- c("chromosome","fork_start","fork_end","readID","read_start","read_end","strand","fork_length_bp","stall_score")
  colnames(rightForks) <- c("chromosome","fork_start","fork_end","readID","read_start","read_end","strand","fork_length_bp","stall_score")
  
  leftForks$fork_direction <- "left"
  rightForks$fork_direction <- "right"
  allForks <- rbind(leftForks, rightForks)
  
  
  #### Origin and Termination ####
  origins <- read.table(paste0(sample_dir,"origins_DNAscent_forkSense.bed"))
  terminations <- read.table(paste0(sample_dir,"terminations_DNAscent_forkSense.bed"))
  
  colnames(origins) <- c("chromosome","origin_start","origin_end","readID","read_start","read_end","strand")
  colnames(terminations) <- c("chromosome","termination_start","termination_end","readID","read_start","read_end","strand")
  
  
  #### Fork Stress Signatures ####
  leftStress <- read.table(paste0(sample_dir,"leftForks_DNAscent_forkSense_stressSignatures.bed"))
  rightStress <- read.table(paste0(sample_dir,"rightForks_DNAscent_forkSense_stressSignatures.bed"))
  
  stress_colnames <- c("chromosome","fork_start","fork_end","readID","read_start","read_end",
                       "strand","fork_length_bp","EdU_length_bp","BrdU_length_bp",
                       "EdU_segment_contamination","EdU_segment_purity",
                       "BrdU_segment_contamination","BrdU_segment_purity",
                       "stall_score")
  
  colnames(leftStress) <- stress_colnames
  colnames(rightStress) <- stress_colnames
  
  leftStress$fork_direction <- "left"
  rightStress$fork_direction <- "right"
  
  allStress <- rbind(leftStress,rightStress)
  
  
  
  #### Save Files ####
  
  write.table(allForks, file=paste0(out_dir,sample,"_allForks.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(allStress, file=paste0(out_dir,sample,"_allStress.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(origins, file=paste0(out_dir,sample,"_origins.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(terminations, file=paste0(out_dir,sample,"_terminations.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(BrdU_segments, file=paste0(out_dir,sample,"_BrdU_segments.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(EdU_segments, file=paste0(out_dir,sample,"_EdU_segments.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  
  
  #### Append samples for collation ####
  
  # Add sample column to every df, in one line
  list2env(lapply(mget(df_names), function(df) { df$sample <- sample; df }), envir = environment())
  
  # Accummulate
  for (nm in df_names) {
    collated[[nm]] <- rbind(collated[[nm]], get(nm))
  }
}

#### Save one combined file per df type, across all samples ####
for (nm in df_names) {
  write.table(collated[[nm]], file = paste0(out_dir, "ALL_samples_", nm, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
