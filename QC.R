# Setting working directory----
setwd("/home/ag-rossi/Timur/")
dir <- file.path("/home/ag-rossi/Timur/")
list.files()

#Loading required packages----
suppressMessages(library(digest))
suppressMessages(library(ShortRead))
suppressMessages(library(knitr))
suppressMessages(library(kableExtra))
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(RColorBrewer))

# Loading functions
source("Static/R/common.R")

#Loading configuration----
resultDir <- file.path("Analysis" , "Results" , "QC")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)
config <- read.table(file = "experiment.txt", header = TRUE, sep = ",")

# Creating the study design, containing all the informations about the samples
studyDesign <- config

studyDesign$readCounts <- file.path("Analysis" , "STAR", 
                                    paste(studyDesign$ID), "Aligned.sortedByCoord.out.bam")
studyDesign$md5 <- lapply(as.character(studyDesign$readCounts), md5sum)

# Check for uniqueness of processed data
sum(duplicated(studyDesign$md5))

# Quality control of aligned reads
flagstatTargets <- file.path("Results" , "Flagstats" ,
                             paste(studyDesign$ID), "flagstats.txt")
file.exists(flagstatTargets)

flagstatRes <- data.frame(matrix(unlist(lapply(flagstatTargets, loadFlagstat)), 
                                 ncol = length(flagstatTargets)), stringsAsFactors = FALSE)

flagstatRes <- as.data.frame(matrix(as.numeric(unlist(flagstatRes)), ncol = length(flagstatTargets)))

mapping_percentage <- t(flagstatRes)
mapping_percentage <- round(as.numeric(mapping_percentage[,8] / mapping_percentage[,2]) * 100, digits = 2 )

flagstatRes <- flagstatRes[-c(6:8),]
flagstatRes <- rbind(flagstatRes, mapping_percentage)

colnames(flagstatRes) <- studyDesign$ID
rownames(flagstatRes) <- c("read mapping" , "primary" , "secondary" , "supplementary" , "duplicates" , "mapping_percentage")

knitr::kable(flagstatRes, caption = "Summary statistics from the STAR alignment", booktabs = TRUE, 
             table.envir = 'table*', linespan = "") %>%
  kableExtra::kable_styling(latex_options = c("hold_position", font_size = 11)) %>%
  add_footnote(c("information from samtools flagstat")) %>%
  scroll_box(width = "2500px", height = "500px")

flagstatRes_t <- as.data.frame(t(flagstatRes))
flagstatRes_t$ID <- row.names(flagstatRes_t)

## Number of reads----
num_read_maps <- ggplot(flagstatRes_t, aes(x = ID ,y = primary))+
  geom_col(fill = "#3b8f51")+
  ggtitle("Number of primary reads per sample")+  
  ylab("Number of reads")+
  scale_y_continuous(limits = c(0, 2.5e8), breaks = seq(0, 2.5e8, by = 1e7))+
  theme(axis.text.x = element_text(angle = 90))
num_read_maps_interactive <- ggplotly(num_read_maps)

htmlwidgets::saveWidget(num_read_maps_interactive, paste0(resultDir, "/Number_primary_reads.html"))

## Percent mapped reads----
percent_mapped_reads <- ggplot(flagstatRes_t, aes(x = ID, y = mapping_percentage))+
  geom_col(fill = "#3b8f51")+
  ggtitle("Percentage of reads mapped")+
  ylab("# mapped reads / # total reads (%)")+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))+
  theme(axis.text.x = element_text(angle = 90))
percent_mapped_reads_interactive <- ggplotly(percent_mapped_reads)

htmlwidgets::saveWidget(percent_mapped_reads_interactive, paste0(resultDir, "/Percent_mapped_reads.html"))
