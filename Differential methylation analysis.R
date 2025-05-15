#loading packages
library(TCGAbiolinks)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(SummarizedExperiment)
library(tidyr)
library(dplyr)
library(sesame)

#creating query
data.met <- GDCquery(project = 'TCGA-LUSC',
                     data.category = 'DNA Methylation',
                     platform = 'Illumina Human Methylation 450',
                     access = 'open',
                     data.type = 'Methylation Beta Value')

#retrieving query
GDCdownload(data.met)
tcga_data <- GDCprepare(data.met)

#inspecting data
dim(tcga_data)
table(tcga_data@colData$vital_status)
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$race)

#reading annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

#removign probes on sex chromosomes to avoid sex-related bias
keep <- !(row.names(tcga_data) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
tcga_data <- tcga_data[keep, ]
rm(keep)

#removing NAs
tcga_data <- subset(tcga_data,subset = (rowSums(is.na(assay(tcga_data))) == 0)) 

#getting clinical data
clinical <- colData(tcga_data)

#reordering
tcga_data <- tcga_data[, rownames(clinical)]
all(rownames(clinical) == colnames(tcga_data))

# replace spaces with "_" in levels of definition column
clinical$definition <-  gsub(" ", "_", clinical$definition)

#creating factor to define levels
clinical$definition <- as.factor(clinical$definition)
clinical$definition <- relevel(clinical$definition, ref = "Solid_Tissue_Normal")


#convert beta values to to mvalues as they are more appropriate for diff methylation analysis
coun <- assay(tcga_data)
mval <- t(apply(coun, 1, function(x) log2(x/(1-x)))) #conversion formula

#analysis
design <- model.matrix(~ definition, data = clinical)

#fitting data to linear model
fit <- lmFit(mval, design)
fit2 <- eBayes(fit)

#saving results
diff_meth = topTable(fit2, coef=ncol(design), sort.by="p",number = nrow(mval), adjust.method = "BH")
write.csv(diff_meth, "DMA.csv")

# Visualization

diff_meth$diffexpressed <- "NO"
diff_meth$diffexpressed[diff_meth$logFC > 2 & diff_meth$adj.P.Val < 0.01] <- "Hyper"
diff_meth$diffexpressed[diff_meth$logFC < -2 & diff_meth$adj.P.Val < 0.01] <- "Hypo"

#volcano plot
ggplot(data = diff_meth, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed, xlab("log2FC"), ylab("p.adj"))) +
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'solid') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'solid') + 
  geom_point(size = 1) +
  scale_color_manual(values = c("#BB2323", "#00009B", "grey" ) ,
                     labels= c("Hyper-Methylated" ,"Hypo-Methylated", "Not significant")) +
  scale_x_continuous(breaks = seq(-10, 10,2)) + scale_y_continuous(breaks = seq(0,350,50))
                     
# annotation of dmgs

annotation <- as.data.frame(ann450k)

merged <- merge(diff_meth, annotation, by = "X") ##X is the column containing probe IDs
