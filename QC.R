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
suppressMessages(library(Rsubread))
suppressMessages(library(AnnotationHub))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))
suppressMessages(library(ashr))
suppressMessages(library(PoiClaClu))
suppressMessages(library(pheatmap))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(VennDiagram))
suppressMessages(library(pathview))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(stringr))

# Loading functions
source("Static/R/common.R")

#Loading configuration----
resultDir <- file.path("Analysis" , "Results" , "QC")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)
config <- read.table(file = "experiment.txt", header = TRUE, sep = ",")

# Creating the study design, containing all the informations about the samples
studyDesign <- config

studyDesign$site_by_age <- paste0(studyDesign$site_of_biopsis, "_" , studyDesign$age_groups)

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

# Assign mapped sequencing reads to specified genomic features----
readCountTarget <- file.path("Analysis" , "STAR",
                             paste(studyDesign$ID), "Aligned.sortedByCoord.out.bam")
readCountTarget_unique <- data.frame(bamFiles = readCountTarget)
readCountTarget_unique$md5sum <- lapply(as.character(readCountTarget_unique$bamFiles), md5sum)


ExternalAnnotation <- file.path("Processing/Annotation/gencode.v32.annotation.gtf")

geneCounts <- featureCounts(files =readCountTarget,
                            annot.ext = ExternalAnnotation,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "exon",
                            GTF.attrType = "gene_id",
                            isLongRead = FALSE,
                            largestOverlap = TRUE,
                            useMetaFeatures = TRUE,
                            nthreads = 16)$counts

colnames(geneCounts) <- studyDesign$ID

# Extraction of the TOP10 genes according to the total number of reads
knitr::kable(geneCounts[order(rowSums(geneCounts), decreasing = TRUE)[1:10],], 
             caption = "Table showing the 10 annotated gene features with the highest number of mapped reads",
             booktabs = TRUE,
             table.envir = 'table*',
             linesep ="") %>%
  kable_styling(latex_options = c("hold_position", font_sire = 11)) %>%
  add_footnote(c("This is raw count data and no normalisation or transformation has been performed")) %>%
  scroll_box(width = "2500px", height = "500px")

#Removal of all genes with zero counts----
geneCounts_nonZeros <- as.data.frame(geneCounts[which(rowSums(geneCounts) > 0), ])
dim(geneCounts)

# Creation of an DEseq2-Data set to store input values
deSeqRaw <- DESeqDataSetFromMatrix(countData = geneCounts_nonZeros, colData = studyDesign, 
                                   design = ~ site_by_age)

# Analysing the variance----
vsd.deSeqRaw <- vst(object = deSeqRaw, blind = TRUE)
topVarGenes <- head(order(rowVars(assay(vsd.deSeqRaw)), decreasing = TRUE), n = 50)
mat <- assay(vsd.deSeqRaw)[topVarGenes,]
mat <- mat - rowMeans(mat)

# Creation of Poisson Distance Matrix----
colors <- colorRampPalette(brewer.pal(4, "Blues"))(255)
poisd <- PoissonDistance(t(geneCounts_nonZeros))
poisd$dd

## Separated by age groups
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(studyDesign$age_groups)
colnames(samplePoisDistMatrix) <- paste(studyDesign$age_groups)

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

## Separated by site
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(studyDesign$site_of_biopsis)
colnames(samplePoisDistMatrix) <- paste(studyDesign$site_of_biopsis)

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# Analyzing differentially expressed genes----
# DE-Analysis based on the binomial distribution
dds <- DESeq(deSeqRaw)
resultsNames(dds)

dds$age_groups <- relevel(dds$age_groups, ref = "young")
dds <- DESeq(dds)
resultsNames(dds)

res_site_of_biopsis <- lfcShrink(dds, coef = 4, type = "apeglm")
plotMA(res_site_of_biopsis,
       colNonSig = "gray32", colSig = "red3",
       log = "x")

res_site_of_biopsis_2 <- results(dds, contrast = c("site_of_biopsis" , "shoulder" , "buttock"))
plotMA(res_site_of_biopsis_2,
       colNonSig = "gray32", colSig = "red3",
       log = "x",
       ylim = c(-3,3))

#Variance analysis of normalized data
vsd <- vst(dds, blind = FALSE) #Estimation of dispersion trend and application of a variance stabilizing transformation

topVarGenesDDS <- head(order(rowVars(assay(vsd)), decreasing = TRUE), n = 25)
mat_dds <- assay(vsd)[topVarGenesDDS,]

mat_dds <- mat_dds - rowMeans(mat_dds)

mat.ens.str <- sub("\\..*","",rownames(mat_dds))
breaksList <- seq(-5,5, by = 0.05)
anno <- as.data.frame(colData(vsd))[, c("site_of_biopsis","age_groups")]
pheatmap(mat_dds, annotation_col = anno, labels_row = mapIds(org.Hs.eg.db, keys = mat.ens.str, column = "SYMBOL",
                                                             keytype = "ENSEMBL", multiVals = "first"),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)

#DEG analysis shoulder vs buttocks sorted by age
#----------young----------#
res_young <- results(dds, contrast = c("site_by_age","shoulder_young","buttock_young"), alpha = 0.1)

DESeq2::plotMA(res_young, ylim =c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("shoulder vs buttocks (young)")

res_young$SYMBOL <- ens.str.SYMBOL(result.table = res_young)
res_young$ENTREZ <- ens.str.ENTREZ(result.table = res_young)

res_young <- res_young[order(abs(res_young$log2FoldChange), decreasing = TRUE) , ]
write.table(x = res_young, file = "Analysis/Results/DE_analysis/res_young.txt" , sep = "\t")

#----------middle----------#
res_middle <- results(dds, contrast = c("site_by_age","shoulder_middle","buttock_middle"), alpha = 0.1)

DESeq2::plotMA(res_middle, ylim =c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("shoulder vs buttocks (middle)")

res_middle$SYMBOL <- ens.str.SYMBOL(result.table = res_middle)
res_middle$ENTREZ <- ens.str.ENTREZ(result.table = res_middle)

res_middle <- res_middle[order(abs(res_middle$log2FoldChange), decreasing = TRUE) , ]
write.table(x = res_middle, file = "Analysis/Results/DE_analysis/res_middle.txt" , sep = "\t")

#----------old----------#
res_old <- results(dds, contrast = c("site_by_age","shoulder_old","buttock_old"), alpha = 0.1)

DESeq2::plotMA(res_old, ylim =c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("shoulder vs buttocks (old)")

res_old$SYMBOL <- ens.str.SYMBOL(result.table = res_old)
res_old$ENTREZ <- ens.str.ENTREZ(result.table = res_old)

res_old <- res_old[order(abs(res_old$log2FoldChange), decreasing = TRUE) , ]
write.table(x = res_old, file = "Analysis/Results/DE_analysis/res_old.txt" , sep = "\t")

# Upregulation and Downregulation - shoulder vs buttocks (age groups)
shoulder_young_total <- data.frame(ENSEMBL_ID = rownames(res_young[abs(res_young$log2FoldChange) > 1.5 , ]))
shoulder_young_up <- data.frame(ENSEMBL_ID = rownames(res_young[res_young$log2FoldChange > 1.5 , ]))
shoulder_young_down <- data.frame(ENSEMBL_ID = rownames(res_young[res_young$log2FoldChange < -1.5 , ]))

shoulder_middle_total <- data.frame(ENSEMBL_ID = rownames(res_middle[abs(res_middle$log2FoldChange) > 1.5 , ]))
shoulder_middle_up <- data.frame(ENSEMBL_ID = rownames(res_middle[res_middle$log2FoldChange > 1.5 , ]))
shoulder_middle_down <- data.frame(ENSEMBL_ID = rownames(res_middle[res_middle$log2FoldChange < -1.5 , ]))

shoulder_old_total <- data.frame(ENSEMBL_ID = rownames(res_old[abs(res_old$log2FoldChange) > 1.5 , ]))
shoulder_old_up <- data.frame(ENSEMBL_ID = rownames(res_old[res_old$log2FoldChange > 1.5 , ]))
shoulder_old_down <- data.frame(ENSEMBL_ID = rownames(res_old[res_old$log2FoldChange < -1.5 , ]))

sheets <- list("shoulder_young_total" = shoulder_young_total,
               "shoulder_young_up" = shoulder_young_up,
               "shoulder_young_down" = shoulder_young_down,
               "shoulder_middle_total" = shoulder_middle_total,
               "shoulder_middle_up" = shoulder_middle_up,
               "shoulder_middle_down" = shoulder_middle_down,
               "shoulder_old_total" = shoulder_old_total,
               "shoulder_old_up" = shoulder_old_up,
               "shoulder_old_down" = shoulder_old_down)

## Creation of Venn diagram
draw.triple.venn(area1 = dim(shoulder_young_total)[1],
                 area2 = dim(shoulder_middle_total)[1],
                 area3 = dim(shoulder_old_total)[1],
                 n12 = length(base::intersect(shoulder_young_total[[1]], shoulder_middle_total[[1]])),
                 n13 = length(base::intersect(shoulder_young_total[[1]], shoulder_old_total[[1]])),
                 n23 = length(base::intersect(shoulder_middle_total[[1]], shoulder_old_total[[1]])),
                 n123 = length(base::intersect(intersect(shoulder_young_total[[1]], shoulder_middle_total[[1]]), shoulder_old_total[[1]])),
                 category = c("young (2436)" , "middle (2407)" , "old (2640)"),
                 euler.d = TRUE,
                 fill = c("#F59B25", "#95E354", "#54A7E3"),
                 alpha = 0.5,
                 fontfamily = rep("Helvetica",7),
                 cat.fontfamily = rep("Helvetica",3),
                 cex = 0.9,
                 cat.cex = 0.9)
writexl::write_xlsx(sheets, "Analysis/Results/DE_analysis/DEG_analysis_L2FC.xlsx")

# Gene Set Enrichment Analysis
#Prepare input
original_gene_list_young <- rownames(res_young)
original_gene_list_young <- res_young[original_gene_list_young, "log2FoldChange"]
names(original_gene_list_young) <-sub("\\..*","",rownames(res_young))
gene_list_young <- na.omit(original_gene_list_young)
gene_list_young <- sort(gene_list_young, decreasing = TRUE)

#Gene Set Enrichment
organism <- org.Hs.eg.db
keytypes(org.Hs.eg.db)

gse_young <- gseGO(geneList = gene_list_young,
             ont = "All",
             keyType = "ENSEMBL",
             minGSSize = 4,
             maxGSSize = 800,
             pvalueCutoff = 0.01,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none",
             seed = TRUE)

dotplot(gse_young, showCategory = 25, x = "GeneRatio", font.size = 12)+
  facet_grid(.~.sign)+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 105))

gse.simple_young <- simplify(x = gseGO(geneList = gene_list_young,
                                       ont = "MF",
                                       keyType = "ENSEMBL",
                                       minGSSize = 4,
                                       maxGSSize = 800,
                                       pvalueCutoff = 0.01,
                                       verbose = TRUE,
                                       OrgDb = organism,
                                       pAdjustMethod = "none",
                                       seed = TRUE),
                             cutoff = 0.7,
                             by = "p.adjust",
                             select_fun = min)

dotplot(gse.simple_young, showCategory = 25, x = "GeneRatio", font.size = 12)+
  facet_grid(.~.sign)+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 105))
