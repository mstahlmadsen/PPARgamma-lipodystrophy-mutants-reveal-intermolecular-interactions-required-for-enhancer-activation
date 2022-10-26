#### ChIP-seq analysis ####
#_________________________________ Read in required packages ________________________________ ####

library(DESeq2)
library(limma)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(gridExtra)
library(ggseqlogo)
library(pROC)

# _________________________________ Analysis HA-ChIPseq ______________________________________####
# Import and annotate data
d.chip <- read.delim("TagCounts_HA_mm10.txt")
size.factors <- read.delim("Size_Factor_HA.txt", header=FALSE)
size.factors_2 <- read.delim("Size_Factor_HA_2.txt", header=FALSE)

rownames(d.chip) <- d.chip$PeakID
colnames(d.chip)[1] <- 'PeakID'
colnames(d.chip)[20:27] <- c('WT_Rep1', 'E379K_Rep1', 'R212Q_Rep1', 'Control_Rep1', 'WT_Rep2', 'E379K_Rep2', 'R212Q_Rep2', 'Control_Rep2')

# Normalize data based on size factor from Homer:
tmp <- d.chip[,20:27]
for(i in 1:dim(tmp)[2]){
  tmp[,i] <- tmp[,i]*size.factors[i,2]
  print(paste0("Using sizefactor for sample ", as.character(size.factors[i,1]), " and multiplying it to sample ", colnames(tmp)[i]))
}
colnames(tmp) <- paste0(colnames(tmp),"_normalized")

# select only peaks with more counts than 35 tags. This is based on visualization in genome browser
tmp <- tmp[ apply(tmp,1,max) > 35,] 
d.chip <- merge(d.chip, tmp, by="row.names")
rownames(d.chip) <- d.chip$PeakID
d.chip$Row.names <- NULL
rm(tmp)

# Remove counts in Pparg-exons as they represents the virus DNA
d.chip <- d.chip[!d.chip$PeakID=="chr6-1167",] 
d.chip <- d.chip[!d.chip$PeakID=="chr6-61",] 
d.chip <- d.chip[!d.chip$PeakID=="chr6-251",] 
d.chip <- d.chip[!d.chip$PeakID=="chr6-334",] 
d.chip <- d.chip[!d.chip$PeakID=="chr6-107",] 
d.chip <- d.chip[!d.chip$PeakID=="chr6-5",] 
d.chip <- d.chip[!d.chip$PeakID=="chr6-34",] 

#### DESeq analysis ####
coldata.chip <- as.data.frame(matrix(nrow=8, ncol=3))
colnames(coldata.chip) <- c("Type", "Replicate", "ID")
coldata.chip[,1] <- c("WT", "WT", "E379K", "E379K", "R212Q", "R212Q", "Control", "Control") 
coldata.chip[,2] <- rep(c("Rep1","Rep2"),4)
coldata.chip[,3] <- paste0(coldata.chip[,1],"_",coldata.chip[,2])
countdata.chip <- d.chip[,c("WT_Rep1", "WT_Rep2", "E379K_Rep1", "E379K_Rep2", "R212Q_Rep1", "R212Q_Rep2")]
size.factors.deseq <- cbind(coldata.chip[,3], 1/size.factors_2$V2)
dds.chip <- DESeqDataSetFromMatrix(countdata.chip, 
                                   coldata.chip[1:6,], 
                                   design = ~ Replicate + Type)
dds.chip <- estimateSizeFactors(dds.chip)
dds.chip$sizeFactor <- as.numeric(size.factors.deseq[1:6,2])
dds.chip <- DESeq(dds.chip)

# PCA plot on top 500 peaks with most variance
p <- as.data.frame(assay(rlog(dds.chip, blind=TRUE)))
p$sd <- apply(as.matrix(p),1,sd)      
p <- p[order(p$sd, decreasing = TRUE),]
p <- p[1:500,] 
p$sd <- NULL

batch <- c("A", "B","A", "B","A", "B")
p.RBE <- removeBatchEffect(p, batch = batch)
p.RBE <- t(p.RBE)
pcadata <- prcomp(p.RBE)
names <-rep(c("WT", "WT","E379K", "E379K","R212Q","R212Q"))
rep <- rep(c("1","2","1","2","1","2"))
ID <- cbind(names, rep, p.RBE)
autoplot(pcadata, data = ID, size=5, shape="rep", scale=0, fill= "names")+
  ggtitle("PCA Top 500 sd peaks") +
  scale_shape_manual(values = c(21:22)) +
  scale_fill_manual(values=c("green4", "blue", "red" )) +
  theme_bw()+ theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1, xlim = c(-20,20), ylim = c(-20,20)) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

# Extracting counts and contrasts
tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "E379K", "WT")))
colnames(tmp) <- paste0(colnames(tmp),"_E379K_WT")
d.chip <- merge(d.chip, tmp, by="row.names")
rownames(d.chip) <- d.chip$PeakID
d.chip$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "WT")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_WT")
d.chip <- merge(d.chip, tmp, by="row.names")
rownames(d.chip) <- d.chip$PeakID
d.chip$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "E379K")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_E379K")
d.chip <- merge(d.chip, tmp, by="row.names")
rownames(d.chip) <- d.chip$PeakID
d.chip$Row.names <- NULL
rm(tmp)

d.chip$WT_average_norm <- (d.chip$WT_Rep1_normalized+d.chip$WT_Rep2_normalized)/2
d.chip$E379K_average_norm <- (d.chip$E379K_Rep1_normalized+d.chip$E379K_Rep2_normalized)/2
d.chip$R212Q_average_norm <- (d.chip$R212Q_Rep1_normalized+d.chip$R212Q_Rep2_normalized)/2
d.chip$Control_average_norm <- (d.chip$Control_Rep1_normalized+d.chip$Control_Rep2_normalized)/2
d.chip$Intensity <- apply(d.chip[,28:35], 1, mean) 

# Generate a bed-file that can be used for counting tags for ATAC-seq, H3K27ac and MED1 ChIPseq within PPARg-peaks that has passed DESeq analysis
write.table(d.chip[,c("Chr", "Start", "End", "PeakID")], col.names = T, row.names = T, sep="\t", quote=F, file="HApeaks_DESeq2.txt") # The format is not bed - converted on server (pos2bed.pl HApeaks_DESeq2.txt > HApeaks_DESeq2.bed). 

#### Output from HA-ChIP DESeq2 analysis ####
write.table(d.chip, col.names = T, row.names = T, sep="\t", quote=F, file="DESeq2_HA_ChIP_TH35.txt")
# cleaning up in R:
rm(list = ls())



# _______________________________ Analysis H3K27ac-ChIPseq ___________________________________####
# Import and annotate data
d.chip.H3K27ac <- read.delim("TagCounts_WTpeaksExt_H3K27ac_mm10.txt")
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
colnames(d.chip.H3K27ac)[1] <- 'PeakID'
colnames(d.chip.H3K27ac)[20:29] <- c('Control_Rep1', 'WT_Rep1', 'E379K_Rep1', 'R212Q_Rep1', 'Control_Rep2', 'WT_Rep2', 'E379K_Rep2', 'R212Q_Rep2', 'Input_Rep1', 'Input_Rep2')
size.factors <- read.delim("Size_Factor_H3K27ac.txt", header=FALSE)
size.factors_2 <- read.delim("Size_Factor_H3K27ac_2.txt", header=FALSE)


# merge with HA-DESeq analysis to filter peaks away that haven't reached the threshold for HA-peaks
d.HA.peaks <- read.delim("HApeaks_DESeq2.bed", header=FALSE)
rownames(d.HA.peaks) <- d.HA.peaks$V4
d.chip.H3K27ac <- merge(d.chip.H3K27ac, d.HA.peaks, by="row.names")
d.chip.H3K27ac$Row.names <- NULL
d.chip.H3K27ac <- d.chip.H3K27ac[1:29]
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
rm(d.HA.peaks)

# Normalize data based on size factor from Homer:
tmp <- d.chip.H3K27ac[,20:29]
for(i in 1:dim(tmp)[2]){
  tmp[,i] <- tmp[,i]*size.factors[i,2]
  print(paste0("Using sizefactor for sample ", as.character(size.factors[i,1]), " and multiplying it to sample ", colnames(tmp)[i]))
}
colnames(tmp) <- paste0(colnames(tmp),"_normalized")

d.chip.H3K27ac <- merge(d.chip.H3K27ac, tmp, by="row.names")
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
d.chip.H3K27ac$Row.names <- NULL
rm(tmp)

#### DESeq2 analysis ####
coldata.chip <- as.data.frame(matrix(nrow=10, ncol=3))
colnames(coldata.chip) <- c("Type", "Replicate", "ID")
coldata.chip[,1] <- c("Control", "Control", "WT", "WT", "E379K", "E379K", "R212Q", "R212Q", "Input", "Input") 
coldata.chip[,2] <- rep(c("Rep1","Rep2"),5)
coldata.chip[,3] <- paste0(coldata.chip[,1],"_",coldata.chip[,2])
countdata.chip <- d.chip.H3K27ac[,c("Control_Rep1", "Control_Rep2", "WT_Rep1", "WT_Rep2", "E379K_Rep1", "E379K_Rep2", "R212Q_Rep1", "R212Q_Rep2")]
size.factors.deseq <- cbind(coldata.chip[,3], 1/size.factors_2$V2)
dds.chip <- DESeqDataSetFromMatrix(countdata.chip, 
                                   coldata.chip[1:8,], 
                                   design = ~ Replicate + Type)
dds.chip <- estimateSizeFactors(dds.chip)
dds.chip$sizeFactor <- as.numeric(size.factors.deseq[1:8,2])
dds.chip <- DESeq(dds.chip)

# PCA plot
pF <- rlog(dds.chip, blind = TRUE)
dataF <- plotPCA(pF, intgroup=c("Type", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(dataF, "percentVar"))
ggplot(dataF, aes(PC1, PC2, color=Type, shape=Replicate)) +
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

# Extracting counts and contrasts
tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "WT", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_WT_Control")
d.chip.H3K27ac <- merge(d.chip.H3K27ac, tmp, by="row.names")
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
d.chip.H3K27ac$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "E379K", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_E379K_Control")
d.chip.H3K27ac <- merge(d.chip.H3K27ac, tmp, by="row.names")
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
d.chip.H3K27ac$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_Control")
d.chip.H3K27ac <- merge(d.chip.H3K27ac, tmp, by="row.names")
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
d.chip.H3K27ac$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "E379K", "WT")))
colnames(tmp) <- paste0(colnames(tmp),"_E379K_WT")
d.chip.H3K27ac <- merge(d.chip.H3K27ac, tmp, by="row.names")
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
d.chip.H3K27ac$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "WT")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_WT")
d.chip.H3K27ac <- merge(d.chip.H3K27ac, tmp, by="row.names")
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
d.chip.H3K27ac$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "E379K")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_E379K")
d.chip.H3K27ac <- merge(d.chip.H3K27ac, tmp, by="row.names")
rownames(d.chip.H3K27ac) <- d.chip.H3K27ac$PeakID
d.chip.H3K27ac$Row.names <- NULL

d.chip.H3K27ac$WT_average_norm <- (d.chip.H3K27ac$WT_Rep1_normalized+d.chip.H3K27ac$WT_Rep2_normalized)/2
d.chip.H3K27ac$E379K_average_norm <- (d.chip.H3K27ac$E379K_Rep1_normalized+d.chip.H3K27ac$E379K_Rep2_normalized)/2
d.chip.H3K27ac$R212Q_average_norm <- (d.chip.H3K27ac$R212Q_Rep1_normalized+d.chip.H3K27ac$R212Q_Rep2_normalized)/2
d.chip.H3K27ac$Control_average_norm <- (d.chip.H3K27ac$Control_Rep1_normalized+d.chip.H3K27ac$Control_Rep2_normalized)/2
d.chip.H3K27ac$Input_average_norm <- (d.chip.H3K27ac$Input_Rep1_normalized+d.chip.H3K27ac$Input_Rep2_normalized)/2
d.chip.H3K27ac$Intensity <- apply(d.chip.H3K27ac[,30:37], 1, mean) 

# Removong 'NA' before further analysis.  
d.chip.H3K27ac$Focus.Ratio.Region.Size <- NULL
#replace NA in columns with padj=NA with 1 (indicating no significant difference)
d.chip.H3K27ac$padj_E379K_WT[is.na(d.chip.H3K27ac$padj_E379K_WT)] <- 1
d.chip.H3K27ac$padj_R212Q_WT[is.na(d.chip.H3K27ac$padj_R212Q_WT)] <- 1
d.chip.H3K27ac$padj_R212Q_E379K[is.na(d.chip.H3K27ac$padj_R212Q_E379K)] <- 1
d.chip.H3K27ac$padj_WT_Control[is.na(d.chip.H3K27ac$padj_WT_Control)] <- 1
d.chip.H3K27ac$padj_E379K_Control[is.na(d.chip.H3K27ac$padj_E379K_Control)] <- 1
d.chip.H3K27ac$padj_R212Q_Control[is.na(d.chip.H3K27ac$padj_R212Q_Control)] <- 1
#replace NA in columns with log2FC=NA with 0 (indicating no significant difference)
d.chip.H3K27ac$log2FoldChange_E379K_WT[is.na(d.chip.H3K27ac$log2FoldChange_E379K_WT)] <- 0
d.chip.H3K27ac$log2FoldChange_R212Q_WT[is.na(d.chip.H3K27ac$log2FoldChange_R212Q_WT)] <- 0
d.chip.H3K27ac$log2FoldChange_R212Q_E379K[is.na(d.chip.H3K27ac$log2FoldChange_R212Q_E379K)] <- 0
d.chip.H3K27ac$log2FoldChange_WT_Control[is.na(d.chip.H3K27ac$log2FoldChange_WT_Control)] <- 0
d.chip.H3K27ac$log2FoldChange_E379K_Control[is.na(d.chip.H3K27ac$log2FoldChange_E379K_Control)] <- 0
d.chip.H3K27ac$log2FoldChange_R212Q_Control[is.na(d.chip.H3K27ac$log2FoldChange_R212Q_Control)] <- 0

#### Output of DESeq analysis ####
write.table(d.chip.H3K27ac, col.names = T, row.names = T, sep="\t", quote=F, file="DESeq2_H3K27_ChIP_HApeaksTH35.txt")
# cleaning up in R:
rm(list = ls())

# ________________________________ Analysis MED1-ChIPseq _____________________________________####
# Import and annotate data
d.chip.MED1 <- read.delim("TagCounts_WTpeaks_MED1_mm10.txt")
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
colnames(d.chip.MED1)[1] <- 'PeakID'
colnames(d.chip.MED1)[20:27] <- c('Control_Rep1', 'WT_Rep1', 'E379K_Rep1', 'R212Q_Rep1', 'Control_Rep2', 'WT_Rep2', 'E379K_Rep2', 'R212Q_Rep2')
size.factors <- read.delim("Size_Factor_MED1.txt", header=FALSE)
size.factors_2 <- read.delim("Size_Factor_MED1_2.txt", header=FALSE)

# merge with HA-DESeq analysis to filter peaks away that haven't reached the threshold for HA-peaks
d.HA.peaks <- read.delim("HApeaks_DESeq2.bed", header=FALSE)
rownames(d.HA.peaks) <- d.HA.peaks$V4
d.chip.MED1 <- merge(d.chip.MED1, d.HA.peaks, by="row.names")
d.chip.MED1$Row.names <- NULL
d.chip.MED1 <- d.chip.MED1[1:27]
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
rm(d.HA.peaks)

# Normalize data based on size factor from Homer:
tmp <- d.chip.MED1[,20:27]
for(i in 1:dim(tmp)[2]){
  tmp[,i] <- tmp[,i]*size.factors[i,2]
  print(paste0("Using sizefactor for sample ", as.character(size.factors[i,1]), " and multiplying it to sample ", colnames(tmp)[i]))
}
colnames(tmp) <- paste0(colnames(tmp),"_normalized")
d.chip.MED1 <- merge(d.chip.MED1, tmp, by="row.names")
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
d.chip.MED1$Row.names <- NULL
rm(tmp)

#### DESeq2 analysis ####
coldata.chip <- as.data.frame(matrix(nrow=8, ncol=3))
colnames(coldata.chip) <- c("Type", "Replicate", "ID")
coldata.chip[,1] <- c("Control", "Control", "WT", "WT", "E379K", "E379K", "R212Q", "R212Q") 
coldata.chip[,2] <- rep(c("Rep1","Rep2"),4)
coldata.chip[,3] <- paste0(coldata.chip[,1],"_",coldata.chip[,2])
countdata.chip <- d.chip.MED1[,c("Control_Rep1", "Control_Rep2", "WT_Rep1", "WT_Rep2", "E379K_Rep1", "E379K_Rep2", "R212Q_Rep1", "R212Q_Rep2")]
size.factors.deseq <- cbind(coldata.chip[,3], 1/size.factors_2$V2)
dds.chip <- DESeqDataSetFromMatrix(countdata.chip, 
                                   coldata.chip, 
                                   design = ~ Replicate + Type)
dds.chip <- estimateSizeFactors(dds.chip)
dds.chip$sizeFactor <- as.numeric(size.factors.deseq[1:8,2])
dds.chip <- DESeq(dds.chip)

# PCA plot on top 500 peaks with most variance
p <- as.data.frame(assay(rlog(dds.chip, blind=TRUE)))
p$sd <- apply(as.matrix(p),1,sd)      
p <- p[order(p$sd, decreasing = TRUE),]
p <- p[1:500,] 
p$sd <- NULL

batch <- c("A", "B","A", "B","A", "B","A", "B")
p.RBE <- removeBatchEffect(p, batch = batch)
p.RBE <- t(p.RBE)
pcadata <- prcomp(p.RBE)
names <-rep(c("Control", "Control", "WT", "WT","E379K", "E379K","R212Q","R212Q"))
rep <- rep(c("1","2","1","2","1","2", "1","2"))
ID <- cbind(names, rep, p.RBE)
autoplot(pcadata, data = ID, size=5, shape="rep", scale=0, fill= "names")+
  ggtitle("PCA Top 500 sd peaks") +
  scale_shape_manual(values = c(21:22)) +
  scale_fill_manual(values=c("grey","green4", "blue", "red" )) +
  theme_bw()+ theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1, xlim = c(-15,15), ylim = c(-15,15)) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

# Extracting counts and contrasts
tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "WT", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_WT_Control")
d.chip.MED1 <- merge(d.chip.MED1, tmp, by="row.names")
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
d.chip.MED1$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "E379K", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_E379K_Control")
d.chip.MED1 <- merge(d.chip.MED1, tmp, by="row.names")
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
d.chip.MED1$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_Control")
d.chip.MED1 <- merge(d.chip.MED1, tmp, by="row.names")
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
d.chip.MED1$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "E379K", "WT")))
colnames(tmp) <- paste0(colnames(tmp),"_E379K_WT")
d.chip.MED1 <- merge(d.chip.MED1, tmp, by="row.names")
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
d.chip.MED1$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "WT")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_WT")
d.chip.MED1 <- merge(d.chip.MED1, tmp, by="row.names")
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
d.chip.MED1$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "E379K")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_E379K")
d.chip.MED1 <- merge(d.chip.MED1, tmp, by="row.names")
rownames(d.chip.MED1) <- d.chip.MED1$PeakID
d.chip.MED1$Row.names <- NULL

d.chip.MED1$WT_average_norm <- (d.chip.MED1$WT_Rep1_normalized+d.chip.MED1$WT_Rep2_normalized)/2
d.chip.MED1$E379K_average_norm <- (d.chip.MED1$E379K_Rep1_normalized+d.chip.MED1$E379K_Rep2_normalized)/2
d.chip.MED1$R212Q_average_norm <- (d.chip.MED1$R212Q_Rep1_normalized+d.chip.MED1$R212Q_Rep2_normalized)/2
d.chip.MED1$Control_average_norm <- (d.chip.MED1$Control_Rep1_normalized+d.chip.MED1$Control_Rep2_normalized)/2
d.chip.MED1$Intensity <- apply(d.chip.MED1[,28:35], 1, mean) 

# Removing 'NA' before further analysis.  
d.chip.MED1$Focus.Ratio.Region.Size <- NULL
#replace NA in columns with padj=NA with 1 (indicating no significant difference)
d.chip.MED1$padj_E379K_WT[is.na(d.chip.MED1$padj_E379K_WT)] <- 1
d.chip.MED1$padj_R212Q_WT[is.na(d.chip.MED1$padj_R212Q_WT)] <- 1
d.chip.MED1$padj_R212Q_E379K[is.na(d.chip.MED1$padj_R212Q_E379K)] <- 1
d.chip.MED1$padj_WT_Control[is.na(d.chip.MED1$padj_WT_Control)] <- 1
d.chip.MED1$padj_E379K_Control[is.na(d.chip.MED1$padj_E379K_Control)] <- 1
d.chip.MED1$padj_R212Q_Control[is.na(d.chip.MED1$padj_R212Q_Control)] <- 1
#replace NA in columns with log2FC=NA with 0 (indicating no significant difference)
d.chip.MED1$log2FoldChange_E379K_WT[is.na(d.chip.MED1$log2FoldChange_E379K_WT)] <- 0
d.chip.MED1$log2FoldChange_R212Q_WT[is.na(d.chip.MED1$log2FoldChange_R212Q_WT)] <- 0
d.chip.MED1$log2FoldChange_R212Q_E379K[is.na(d.chip.MED1$log2FoldChange_R212Q_E379K)] <- 0
d.chip.MED1$log2FoldChange_WT_Control[is.na(d.chip.MED1$log2FoldChange_WT_Control)] <- 0
d.chip.MED1$log2FoldChange_E379K_Control[is.na(d.chip.MED1$log2FoldChange_E379K_Control)] <- 0
d.chip.MED1$log2FoldChange_R212Q_Control[is.na(d.chip.MED1$log2FoldChange_R212Q_Control)] <- 0

#### Output from Med1-ChIP DESeq2 analysis ####
write.table(d.chip.MED1, col.names = T, row.names = T, sep="\t", quote=F, file="DESeq2_MED1_ChIP_HApeaksTH35.txt")
# cleaning up in R:
rm(list = ls())

# ___________________________________ Analysis ATACseq _______________________________________####
# Import and annotate data
d.chip.ATAC <- read.delim("TagCounts_WTpeaks_ATAC_mm10.txt")
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
colnames(d.chip.ATAC)[1] <- 'PeakID'
colnames(d.chip.ATAC)[20:27] <- c('Control_Rep1', 'Control_Rep2','E379K_Rep1', 'E379K_Rep2', 'R212Q_Rep1', 'R212Q_Rep2', 'WT_Rep1', 'WT_Rep2')
size.factors <- read.delim("Size_Factor_ATAC.txt", header=FALSE)
size.factors_2 <- read.delim("Size_Factor_ATAC_2.txt", header=FALSE)

# merge with HA-DESeq analysis to filter peaks away that haven't reached the threshold for HA-peaks
d.HA.peaks <- read.delim("HApeaks_DESeq2.bed", header=FALSE)
rownames(d.HA.peaks) <- d.HA.peaks$V4
d.chip.ATAC <- merge(d.chip.ATAC, d.HA.peaks, by="row.names")
d.chip.ATAC$Row.names <- NULL
d.chip.ATAC <- d.chip.ATAC[1:27]
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
rm(d.HA.peaks)

# Normalize data based on size factor from Homer:
tmp <- d.chip.ATAC[,20:27]
for(i in 1:dim(tmp)[2]){
  tmp[,i] <- tmp[,i]*size.factors[i,2]
  print(paste0("Using sizefactor for sample ", as.character(size.factors[i,1]), " and multiplying it to sample ", colnames(tmp)[i]))
}
colnames(tmp) <- paste0(colnames(tmp),"_normalized")

d.chip.ATAC <- merge(d.chip.ATAC, tmp, by="row.names")
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
d.chip.ATAC$Row.names <- NULL
rm(tmp)

#### DESeq2 analysis ####
coldata.chip <- as.data.frame(matrix(nrow=8, ncol=3))
colnames(coldata.chip) <- c("Type", "Replicate", "ID")
coldata.chip[,1] <- c("Control", "Control", "WT", "WT", "E379K", "E379K", "R212Q", "R212Q") 
coldata.chip[,2] <- rep(c("Rep1","Rep2"),4)
coldata.chip[,3] <- paste0(coldata.chip[,1],"_",coldata.chip[,2])
countdata.chip <- d.chip.ATAC[,c("Control_Rep1", "Control_Rep2", "WT_Rep1", "WT_Rep2", "E379K_Rep1", "E379K_Rep2", "R212Q_Rep1", "R212Q_Rep2")]
size.factors.deseq <- cbind(coldata.chip[,3], 1/size.factors_2$V2)
dds.chip <- DESeqDataSetFromMatrix(countdata.chip, 
                                   coldata.chip, 
                                   design = ~ Replicate + Type)
dds.chip <- estimateSizeFactors(dds.chip)
dds.chip$sizeFactor <- as.numeric(size.factors.deseq[1:8,2])
dds.chip <- DESeq(dds.chip)

# PCA plot on top 500 peaks with most variance
p <- as.data.frame(assay(rlog(dds.chip, blind=TRUE)))
p$sd <- apply(as.matrix(p),1,sd)      
p <- p[order(p$sd, decreasing = TRUE),]
p <- p[1:500,] 
p$sd <- NULL

batch <- c("A", "B","A", "B","A", "B","A", "B")
p.RBE <- removeBatchEffect(p, batch = batch)
p.RBE <- t(p.RBE)
pcadata <- prcomp(p.RBE)
names <-rep(c("Control", "Control", "WT", "WT","E379K", "E379K","R212Q","R212Q"))
rep <- rep(c("1","2","1","2","1","2", "1","2"))
ID <- cbind(names, rep, p.RBE)
autoplot(pcadata, data = ID, size=5, shape="rep", scale=0, fill= "names")+
  ggtitle("PCA Top 500 sd peaks") +
  scale_shape_manual(values = c(21:22)) +
  scale_fill_manual(values=c("grey","green4", "blue", "red" )) +
  theme_bw()+ theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1, xlim = c(-7,7), ylim = c(-7,7)) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

# Extracting counts and contrasts
tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "WT", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_WT_Control")
d.chip.ATAC <- merge(d.chip.ATAC, tmp, by="row.names")
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
d.chip.ATAC$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "E379K", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_E379K_Control")
d.chip.ATAC <- merge(d.chip.ATAC, tmp, by="row.names")
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
d.chip.ATAC$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "Control")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_Control")
d.chip.ATAC <- merge(d.chip.ATAC, tmp, by="row.names")
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
d.chip.ATAC$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "E379K", "WT")))
colnames(tmp) <- paste0(colnames(tmp),"_E379K_WT")
d.chip.ATAC <- merge(d.chip.ATAC, tmp, by="row.names")
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
d.chip.ATAC$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "WT")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_WT")
d.chip.ATAC <- merge(d.chip.ATAC, tmp, by="row.names")
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
d.chip.ATAC$Row.names <- NULL

tmp <- as.data.frame(results(dds.chip, contrast = c("Type", "R212Q", "E379K")))
colnames(tmp) <- paste0(colnames(tmp),"_R212Q_E379K")
d.chip.ATAC <- merge(d.chip.ATAC, tmp, by="row.names")
rownames(d.chip.ATAC) <- d.chip.ATAC$PeakID
d.chip.ATAC$Row.names <- NULL

d.chip.ATAC$WT_average_norm <- (d.chip.ATAC$WT_Rep1_normalized+d.chip.ATAC$WT_Rep2_normalized)/2
d.chip.ATAC$E379K_average_norm <- (d.chip.ATAC$E379K_Rep1_normalized+d.chip.ATAC$E379K_Rep2_normalized)/2
d.chip.ATAC$R212Q_average_norm <- (d.chip.ATAC$R212Q_Rep1_normalized+d.chip.ATAC$R212Q_Rep2_normalized)/2
d.chip.ATAC$Control_average_norm <- (d.chip.ATAC$Control_Rep1_normalized+d.chip.ATAC$Control_Rep2_normalized)/2
d.chip.ATAC$Intensity <- apply(d.chip.ATAC[,28:35], 1, mean) 

# Removing 'NA' before further analysis.  
d.chip.ATAC$Focus.Ratio.Region.Size <- NULL
#replace NA in columns with padj=NA with 1 (indicating no significant difference)
d.chip.ATAC$padj_E379K_WT[is.na(d.chip.ATAC$padj_E379K_WT)] <- 1
d.chip.ATAC$padj_R212Q_WT[is.na(d.chip.ATAC$padj_R212Q_WT)] <- 1
d.chip.ATAC$padj_R212Q_E379K[is.na(d.chip.ATAC$padj_R212Q_E379K)] <- 1
d.chip.ATAC$padj_WT_Control[is.na(d.chip.ATAC$padj_WT_Control)] <- 1
d.chip.ATAC$padj_E379K_Control[is.na(d.chip.ATAC$padj_E379K_Control)] <- 1
d.chip.ATAC$padj_R212Q_Control[is.na(d.chip.ATAC$padj_R212Q_Control)] <- 1
#replace NA in columns with log2FC=NA with 0 (indicating no significant difference)
d.chip.ATAC$log2FoldChange_E379K_WT[is.na(d.chip.ATAC$log2FoldChange_E379K_WT)] <- 0
d.chip.ATAC$log2FoldChange_R212Q_WT[is.na(d.chip.ATAC$log2FoldChange_R212Q_WT)] <- 0
d.chip.ATAC$log2FoldChange_R212Q_E379K[is.na(d.chip.ATAC$log2FoldChange_R212Q_E379K)] <- 0
d.chip.ATAC$log2FoldChange_WT_Control[is.na(d.chip.ATAC$log2FoldChange_WT_Control)] <- 0
d.chip.ATAC$log2FoldChange_E379K_Control[is.na(d.chip.ATAC$log2FoldChange_E379K_Control)] <- 0
d.chip.ATAC$log2FoldChange_R212Q_Control[is.na(d.chip.ATAC$log2FoldChange_R212Q_Control)] <- 0

#### Output of DESeq analysis ####
write.table(d.chip.ATAC, col.names = T, row.names = T, sep="\t", quote=F, file="DESeq2_ATAC_ChIP_HApeaksTH35.txt")
# cleaning up in R:
rm(list = ls())

# _______________________________ Motif analysis - known motif _______________________________####
# Annotate JASPAR PPRE to all PPARg-binding sites.
# Run on server:    annotatePeaks.pl HA_peaks_35.txt mm10 -noann -nogene -size 200 -mscore -m MA0065.2_7.motif > HA_peaks_35_200bp_motifscore.txt
# Import and annotate data
Motif <- read.delim("HA_peaks_35_200bp_motifscore.txt")
colnames(Motif)[1] <- 'PeakID'
rownames(Motif) <- Motif$PeakID
Motif_2 <- as.data.frame(Motif[,10])
row.names(Motif_2) <- row.names(Motif)
colnames(Motif_2) <- c("Jaspar_PPRE")

# _______________________________ Analysis PPARg-CEBPa synergy _______________________________####
# Data from DOI: 10.1128/MCB.01344-13
d.syn <- read.delim("C:/Users/msm/OneDrive - Syddansk Universitet/Projects/PPARg-CEBPa_synergy/Reanalysis_2020_forFPLD3story/TagCount_PPARgCEBPa_mm10_PPARgWTTH35peaks.txt")
rownames(d.syn) <- d.syn$PeakID
colnames(d.syn)[1] <- 'PeakID'
colnames(d.syn)[20:26] <- c('SC14AA_CEBPa', 'SC14AA_Control', 'SC14AA_PPARgCEBPa', 'input', 'H100_Control', 'H100_PPARgCEBPa', 'H100_PPARg')
d.syn <- d.syn[,c(1:19, 23, 21,20, 22, 24, 26, 25)]

# Determine if PPARg and C/EBPa from 'old' data are binding to HA-PPARg-peaks (same criteria as in MCB paper):
tmp <- d.syn[,c(20,24:26)]
tmp <- tmp[ apply(tmp,1,max) >= 20,] #filter away peaks with less than 20 tags
tmp <- tmp[tmp$H100_PPARg/tmp$input >=4 |tmp$H100_PPARgCEBPa/tmp$input>=4,] # filter away sites where binding is less than 4-fold the input-signal
PPARg_bind <- merge(d.syn, tmp, by="row.names")
rownames(PPARg_bind) <- PPARg_bind$PeakID
PPARg_bind$Row.names <- NULL
PPARg_bind <- PPARg_bind[,1:26]
colnames(PPARg_bind)[20] <- 'input'
colnames(PPARg_bind)[24:26] <- c('H100_Control', 'H100_PPARg', 'H100_PPARgCEBPa')
rm(tmp)

tmp <- d.syn[,c(20:23)]
tmp <- tmp[ apply(tmp,1,max) >= 20,] #filter away peaks with less than 20 tags
tmp <- tmp[tmp$SC14AA_CEBPa/tmp$input >=4 |tmp$SC14AA_PPARgCEBPa/tmp$input>=4,] # filter away sites where binding is less than 4-fold the input-signal
CEBPa_bind <- merge(d.syn, tmp, by="row.names")
rownames(CEBPa_bind) <- CEBPa_bind$PeakID
CEBPa_bind$Row.names <- NULL
CEBPa_bind <- CEBPa_bind[,1:26]
colnames(CEBPa_bind)[20] <- 'input'
colnames(CEBPa_bind)[21:23] <- c('SC14AA_Control', 'SC14AA_CEBPa', 'SC14AA_PPARgCEBPa')
rm(tmp)

d.syn$PPARg_bind <- 0
d.syn$PPARg_bind <- ifelse(d.syn$PeakID %in% PPARg_bind$PeakID, 1,0)
d.syn$CEBPa_bind <- 0
d.syn$CEBPa_bind <- ifelse(d.syn$PeakID %in% CEBPa_bind$PeakID, 1,0)

# Defining enhancers that gain, lose or have constitutive binding dependent on TF-coexpression (same criteria as in MCB paper)
PPARg_gain <- d.syn[d.syn$PPARg_bind ==1 & d.syn$H100_PPARgCEBPa/d.syn$H100_PPARg>=2.5,]
PPARg_lost <- d.syn[d.syn$PPARg_bind ==1 & d.syn$H100_PPARgCEBPa/d.syn$H100_PPARg<=0.4,]
PPARg_const <- d.syn[d.syn$PPARg_bind ==1 & d.syn$H100_PPARgCEBPa/d.syn$H100_PPARg>0.4 & d.syn$H100_PPARgCEBPa/d.syn$H100_PPARg<2.5,]

CEBPa_gain <- d.syn[d.syn$CEBPa_bind ==1 & d.syn$SC14AA_PPARgCEBPa/d.syn$SC14AA_CEBPa>=2.5,]
CEBPa_lost <- d.syn[d.syn$CEBPa_bind ==1 & d.syn$SC14AA_PPARgCEBPa/d.syn$SC14AA_CEBPa<=0.4,]
CEBPa_const <- d.syn[d.syn$CEBPa_bind ==1 & d.syn$SC14AA_PPARgCEBPa/d.syn$SC14AA_CEBPa>0.4 & d.syn$SC14AA_PPARgCEBPa/d.syn$SC14AA_CEBPa<2.5,]

# Merging synergy-groups into master_df
# PPARg_eff = effect of co-introducing PPARg
# CEBPa_eff = effect of co-introducing CEBPa

PPARg_const$CEBPa_eff <- 'Const'
PPARg_gain$CEBPa_eff <- 'Gain'
PPARg_lost$CEBPa_eff <- 'Lost'
No_PPARg <- d.syn[!(d.syn$PeakID %in% PPARg_bind$PeakID),]
No_PPARg$CEBPa_eff <- '0'
PPARg_bind <- rbind(PPARg_const, PPARg_gain, PPARg_lost, No_PPARg)

CEBPa_const$PPARg_eff <- 'Const'
CEBPa_gain$PPARg_eff <- 'Gain'
CEBPa_lost$PPARg_eff <- 'Lost'
No_CEBPa <- d.syn[!(d.syn$PeakID %in% CEBPa_bind$PeakID),]
No_CEBPa$PPARg_eff <- '0'
CEBPa_bind <- rbind(CEBPa_const, CEBPa_gain, CEBPa_lost, No_CEBPa)

tmp <- PPARg_bind[,c(1,29)]  
d.syn <- merge(d.syn, tmp, by='PeakID')
tmp <- CEBPa_bind[,c(1,29)]
d.syn <- merge(d.syn, tmp, by='PeakID')
rownames(d.syn) <- d.syn$PeakID
rm(tmp)

# Taking required information:
d.syn <- d.syn[,c(1, 20:30)]
write.table(d.syn, col.names = T, row.names = T, sep="\t", quote=F, file="HA_peaks_35_PPARgCEBPa_syn.txt")

# _____________________________ Combining data in one data frame _____________________________####
# Importing data
HA_data <- read.delim("DESeq2_HA_ChIP_TH35.txt")
MED1_data <- read.delim("DESeq2_MED1_ChIP_HApeaksTH35.txt")
H3K27ac_data <- read.delim("DESeq2_H3K27_ChIP_HApeaksTH35.txt")
ATAC_data <- read.delim("DESeq2_ATAC_ChIP_HApeaksTH35.txt")

# Combining data frames
tmp <- MED1_data[,c(27:75)]
colnames(tmp) <- paste0("MED1_", colnames(tmp))
d.chip.master <- merge(HA_data, tmp, by="row.names")
rownames(d.chip.master) <- d.chip.master$PeakID
d.chip.master$Row.names <- NULL
rm(MED1_data)
rm(HA_data)

tmp <- H3K27ac_data[,c(29:80)]
colnames(tmp) <- paste0("H3K27ac_", colnames(tmp))
d.chip.master <- merge(d.chip.master, tmp, by="row.names")
rownames(d.chip.master) <- d.chip.master$PeakID
d.chip.master$Row.names <- NULL
rm(H3K27ac_data)

tmp <- ATAC_data[,c(27:75)]
colnames(tmp) <- paste0("ATAC_", colnames(tmp))
d.chip.master <- merge(d.chip.master, tmp, by="row.names")
rownames(d.chip.master) <- d.chip.master$PeakID
d.chip.master$Row.names <- NULL
rm(ATAC_data)
rm(tmp)

d.chip.master <- merge(d.chip.master, Motif_2, by="row.names")
rownames(d.chip.master) <- d.chip.master$PeakID
d.chip.master$Row.names <- NULL
rm(Motif, Motif_2)

# _______________________________________  Figure 6 __________________________________________####
#### Defining enhancer groups ####
FDR <- 0.1
L2FC <- 0
L2FC2 <- log2(1.25)

# H3K27ac responsive enhancers
WT.gain.H3K27ac <- d.chip.master[d.chip.master$H3K27ac_padj_WT_Control < FDR & d.chip.master$H3K27ac_log2FoldChange_WT_Control > L2FC, ]
WT.lose.H3K27ac <- d.chip.master[d.chip.master$H3K27ac_padj_WT_Control < FDR & d.chip.master$H3K27ac_log2FoldChange_WT_Control < L2FC, ]
WT.constant.H3K27ac <- d.chip.master[!(d.chip.master$PeakID %in% WT.gain.H3K27ac$PeakID) & !(d.chip.master$PeakID %in% WT.lose.H3K27ac$PeakID),]

# MED1 responsive enhancers
WT.gain.MED1 <- d.chip.master[d.chip.master$MED1_padj_WT_Control < FDR & d.chip.master$MED1_log2FoldChange_WT_Control > L2FC, ]
WT.lose.MED1 <- d.chip.master[d.chip.master$MED1_padj_WT_Control < FDR & d.chip.master$MED1_log2FoldChange_WT_Control < L2FC, ]
WT.constant.MED1 <- d.chip.master[!(d.chip.master$PeakID %in% WT.gain.MED1$PeakID) & !(d.chip.master$PeakID %in% WT.lose.MED1$PeakID),]

# H3K27ac and MED1 responsonsive enhancers
WT.gain.H3K27ac.MED1 <- d.chip.master[(d.chip.master$PeakID %in% WT.gain.H3K27ac$PeakID) & (d.chip.master$PeakID %in% WT.gain.MED1$PeakID),]
WT.gain.H3K27acOrMED1 <- d.chip.master[(d.chip.master$PeakID %in% WT.gain.H3K27ac$PeakID) | (d.chip.master$PeakID %in% WT.gain.MED1$PeakID),]
target.enhancers <- WT.gain.H3K27acOrMED1
rm(WT.gain.H3K27acOrMED1)

# defining enhancers that does not gain H3K27ac or MED1 signal:
non_activated.enhancers <- d.chip.master[(!d.chip.master$PeakID %in% target.enhancers$PeakID) & !(d.chip.master$PeakID %in% WT.lose.H3K27ac$PeakID) & !(d.chip.master$PeakID %in% WT.lose.MED1$PeakID),]

# enhancers that does not change H3K27ac or MED1 (more strict definition):
constant.enhancers <- non_activated.enhancers[non_activated.enhancers$H3K27ac_log2FoldChange_WT_Control <L2FC2 
                                              & non_activated.enhancers$H3K27ac_log2FoldChange_WT_Control > -(L2FC2)
                                              & non_activated.enhancers$MED1_log2FoldChange_WT_Control <L2FC2 
                                              & non_activated.enhancers$MED1_log2FoldChange_WT_Control > -(L2FC2),] 

#### Fig 6.b ####
a <- ggplot(WT.constant.H3K27ac, aes(x=log(H3K27ac_WT_average_norm,10), y=H3K27ac_log2FoldChange_WT_Control))+
        geom_point(size =3, alpha=0.6, shape=16, colour="grey")+
        ggtitle("H3K27ac response to PPARg2WT expression") + 
        ylim(-3.5,3.5)+
        geom_point(data=WT.gain.H3K27ac, aes(x=log(H3K27ac_WT_average_norm,10), y=H3K27ac_log2FoldChange_WT_Control), colour="red3", size =3, alpha=0.6, shape=16)+
        geom_point(data=WT.lose.H3K27ac, aes(x=log(H3K27ac_WT_average_norm,10), y=H3K27ac_log2FoldChange_WT_Control), colour="black", size =3, alpha=0.6, shape=16)+
        #geom_point(data=constant.enhancers, aes(x=log(H3K27ac_WT_average_norm,10), y=H3K27ac_log2FoldChange_WT_Control), colour="blue", size =3, alpha=0.7, shape=1)+
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
        geom_hline(yintercept=0)

b <- ggplot(WT.constant.MED1, aes(x=log(MED1_WT_average_norm,10), y=MED1_log2FoldChange_WT_Control))+
        geom_point(size =3, alpha=0.6, shape=16, colour="grey")+
        ggtitle("MED1 response to PPARg2WT expression") + 
        ylim(-8,8)+
        geom_point(data=WT.gain.MED1, aes(x=log(MED1_WT_average_norm,10), y=MED1_log2FoldChange_WT_Control), colour="red", size =3, alpha=0.6, shape=16)+
        geom_point(data=WT.lose.MED1, aes(x=log(MED1_WT_average_norm,10), y=MED1_log2FoldChange_WT_Control), colour="black", size =3, alpha=0.6, shape=16)+
        #geom_point(data=constant.enhancers, aes(x=log(MED1_WT_average_norm,10), y=MED1_log2FoldChange_WT_Control), colour="blue", size =3, alpha=0.7, shape=1)+
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
        geom_hline(yintercept=0)

pdf('Fig6b_HA_PPARg_targetEnhancers_MAplot.pdf', width=10, height = 5)
grid.arrange(a, b, ncol = 2, nrow = 1)
dev.off()
rm(a,b)

#### Fig 6.c ####

print(paste0(nrow(WT.gain.MED1), " peaks gain MED1 occupancy upon PPARg expression"))
print(paste0(nrow(WT.gain.H3K27ac), " peaks gain H3K27ac upon PPARg expression"))
print(paste0(nrow(WT.gain.H3K27ac.MED1), " peaks gain both H3K27ac and MED1 occupancy upon PPARg expression"))
print(paste0(nrow(target.enhancers), " peaks are PPARg-target enhancers"))
# Venn diagram plotted using http://eulerr.co/

#### Fig 6.d ####
pdf('Fig6D_target-enhancers_boxplot_HA-ChIP.pdf', width=5, height = 5)
par(pty="s")
boxplot(non_activated.enhancers$WT_average_norm,
        target.enhancers$WT_average_norm,
        outline=F,
        notch=TRUE,
        ylim=c(0, 600),
        main= "PPARg2-WT\nChIP-seq signal",
        ylab='PPARg2-WT normalized counts',
        names=c("non-activated", "activated"),
        las=2,
        col=c("grey", "red"))
dev.off()
wilcox.test(non_activated.enhancers$WT_average_norm, target.enhancers$WT_average_norm)$p.value
print(paste0(nrow(non_activated.enhancers), " enhancers are not activated by PPARg"))
print(paste0(nrow(target.enhancers), " enhancers are activated by PPARg"))


#### Fig 6.e ####
pdf('Fig6E_target-enhancers_boxplot_Jaspar.pdf', width=5, height = 5)
boxplot(non_activated.enhancers$Jaspar_PPRE,
        target.enhancers$Jaspar_PPRE,
        outline=F,
        notch=TRUE,
        ylim=c(-1, 16),
        main= "JASPAR PPRE\nmotif strength",
        ylab='Log odds motif score',
        names=c("non-activated",  "activated"),
        las=2,
        col=c("grey","red"))
dev.off()
wilcox.test(non_activated.enhancers$Jaspar_PPRE, target.enhancers$Jaspar_PPRE)$p.value

#### De novo motif analysis ####
# write position file for PPARg-target enhancers:
write.table(target.enhancers[,c("PeakID", "Chr", "Start", "End", "Strand")], col.names = T, row.names = F, sep="\t", quote=F, file="target_enhancers.txt") 
# Perform doNovo motif search within PPARg-target enhancers:
# Run on server:    findMotifsGenome.pl target_enhancers.txt mm10 Target_enhancers_Denovo.motifs/ -len 15,16,17 -size 200 -S 15 -bits
# The top scoring motif from this analysis resembles the PPRE in the revers orientation. Annotate this motif to all PPARg-peaks passing DESeq2 analysis. 
# Run on server:    annotatePeaks.pl HA_peaks_35.txt mm10 -noann -nogene -size 200 -mscore -m motif1RV.motif > HA_peaks_35_deNovoMotifscore_TE_200bp.txt
# Importing and annotating data
DeNovo.Motif.TE <- read.delim("HA_peaks_35_deNovoMotifscore_TE_200bp.txt")
colnames(DeNovo.Motif.TE)[1] <- 'PeakID'
rownames(DeNovo.Motif.TE) <- DeNovo.Motif.TE$PeakID
DeNovo.Motif.TE_2 <- as.data.frame(DeNovo.Motif.TE[,10])
row.names(DeNovo.Motif.TE_2) <- row.names(DeNovo.Motif.TE)
colnames(DeNovo.Motif.TE_2) <- c("DeNovo_PPRE_TE")

# Integrating into the dataframe 'target.enhancers'
target.enhancers <- merge(target.enhancers, DeNovo.Motif.TE_2, by="row.names")
rownames(target.enhancers) <- target.enhancers$PeakID
target.enhancers$Row.names <- NULL
rm(DeNovo.Motif.TE, DeNovo.Motif.TE_2)

#### Fig 6.f ####
#JASPAR PPRE
PWMforward <- read.delim("MA0065.2_7.motif", skip=1, header=F)
Logo <- t(PWMforward)
row.names(Logo) <- c("A","C", "G", "T")
pdf('JASPAR_motif_PWM.pdf', width=10, height = 5)
ggplot() + geom_logo(Logo) + theme_logo()
dev.off()

# DeNovo motif
PWMforward <- read.delim("motif1RV.motif", skip=1, header=F)
Logo <- t(PWMforward)
row.names(Logo) <- c("A","C", "G", "T")
pdf('Motifs_PWM.pdf', width=10, height = 5)
ggplot() + geom_logo(Logo) + theme_logo()
dev.off()

#### Fig 6.g ####
pdf('Fig5G_mut_sensitivety_target-enhancers_boxplot.pdf', width=10, height = 15)
par(mfrow=c(3,2), pty="s")
boxplot(non_activated.enhancers$log2FoldChange_E379K_WT,
        target.enhancers$log2FoldChange_E379K_WT,
        outline=F,
        notch=TRUE,
        main= "E379K",
        ylab='PPARgChIP-seq \nlog2FoldChange(mutant vs. WT)',
        names=c("non-activated",  "activated"),
        ylim=c(-4.5, 3),
        las=2,
        col=c("grey",  "red"))
abline(h=0, lty=2)

boxplot(non_activated.enhancers$log2FoldChange_R212Q_WT,
        target.enhancers$log2FoldChange_R212Q_WT,
        outline=F,
        notch=TRUE,
        main= "R212Q",
        ylab='PPARgChIP-seq \nlog2FoldChange(mutant vs. WT)',
        names=c("non-activated",  "activated"),
        ylim=c(-4.5, 3),
        las=2,
        col=c("grey", "red"))
abline(h=0, lty=2)

boxplot(non_activated.enhancers$H3K27ac_log2FoldChange_E379K_WT,
        target.enhancers$H3K27ac_log2FoldChange_E379K_WT,
        outline=F,
        notch=TRUE,
        main= "E379K",
        ylab='H3K27ac ChIP-seq \nlog2FoldChange(mutant vs. WT)',
        names=c("non-activated",  "activated"),
        ylim=c(-2, 1),
        las=2,
        col=c("grey",  "red"))
abline(h=0, lty=2)

boxplot(non_activated.enhancers$H3K27ac_log2FoldChange_R212Q_WT,
        target.enhancers$H3K27ac_log2FoldChange_R212Q_WT,
        outline=F,
        notch=TRUE,
        main= "R212Q",
        ylab='H3K27ac ChIP-seq \nlog2FoldChange(mutant vs. WT)',
        names=c("non-activated", "activated"),
        ylim=c(-2, 1),
        las=2,
        col=c("grey", "red"))
abline(h=0, lty=2)

boxplot(non_activated.enhancers$MED1_log2FoldChange_E379K_WT,
        target.enhancers$MED1_log2FoldChange_E379K_WT,
        outline=F,
        notch=TRUE,
        main= "E379K",
        ylab='MED1 ChIP-seq \nlog2FoldChange(mutant vs. WT)',
        names=c("non-activated", "activated"),
        ylim=c(-3, 1.5),
        las=2,
        col=c("grey", "red"))
abline(h=0, lty=2)

boxplot(non_activated.enhancers$MED1_log2FoldChange_R212Q_WT,
        target.enhancers$MED1_log2FoldChange_R212Q_WT,
        outline=F,
        notch=TRUE,
        main= "R212Q",
        ylab='MED1 ChIP-seq \nlog2FoldChange(mutant vs. WT)',
        names=c("non-activated", "activated"),
        ylim=c(-3, 1.5),
        las=2,
        col=c("grey", "red"))
abline(h=0, lty=2)
dev.off()

wilcox.test(non_activated.enhancers$log2FoldChange_E379K_WT, target.enhancers$log2FoldChange_E379K_WT)$p.value
wilcox.test(non_activated.enhancers$log2FoldChange_R212Q_WT, target.enhancers$log2FoldChange_R212Q_WT)$p.value

wilcox.test(non_activated.enhancers$H3K27ac_log2FoldChange_E379K_WT, target.enhancers$H3K27ac_log2FoldChange_E379K_WT)$p.value
wilcox.test(non_activated.enhancers$H3K27ac_log2FoldChange_R212Q_WT, target.enhancers$H3K27ac_log2FoldChange_R212Q_WT)$p.value

wilcox.test(non_activated.enhancers$MED1_log2FoldChange_E379K_WT, target.enhancers$MED1_log2FoldChange_E379K_WT)$p.value
wilcox.test(non_activated.enhancers$MED1_log2FoldChange_R212Q_WT, target.enhancers$MED1_log2FoldChange_R212Q_WT)$p.value












# _______________________________________  Figure 7 __________________________________________####
#### Defining mutation-sensitive enhancers ####
target.enhancers$Position <- paste0(target.enhancers[,16],"_",target.enhancers[,10])

FDR2<- 0.05
L2FC3<-log2(1.25)

E379K_induced <- target.enhancers[ target.enhancers$padj_E379K_WT < FDR2 & target.enhancers$log2FoldChange_E379K_WT >L2FC3 ,]
E379K_repressed <- target.enhancers[target.enhancers$padj_E379K_WT < FDR2 & target.enhancers$log2FoldChange_E379K_WT <L2FC3,]
E379K_constant <- target.enhancers[!(target.enhancers$PeakID %in% E379K_induced$PeakID) & !(target.enhancers$PeakID %in% E379K_repressed$PeakID),] # loose definition, all target-enhancers not significantly induced or repressed
non.changing.E379K <- target.enhancers[target.enhancers$padj_E379K_WT >0.1 & abs(target.enhancers$log2FoldChange_E379K_WT) < 0.25, ] # more stringent, peaks changing less than 2-fold

R212Q_induced <- target.enhancers[ target.enhancers$padj_R212Q_WT < FDR2 & target.enhancers$log2FoldChange_R212Q_WT >L2FC3 ,]
R212Q_repressed <- target.enhancers[target.enhancers$padj_R212Q_WT < FDR2 & target.enhancers$log2FoldChange_R212Q_WT <L2FC3,]
R212Q_constant <- target.enhancers[!(target.enhancers$PeakID %in% R212Q_induced$PeakID) & !(target.enhancers$PeakID %in% R212Q_repressed$PeakID),]
non.changing.R212Q <- target.enhancers[target.enhancers$padj_R212Q_WT >0.1 & abs(target.enhancers$log2FoldChange_R212Q_WT) < 0.25, ]


mut.rep <- target.enhancers[ target.enhancers$padj_E379K_WT < FDR2 & target.enhancers$log2FoldChange_E379K_WT < L2FC3 |
                               target.enhancers$padj_R212Q_WT < FDR2 & target.enhancers$log2FoldChange_R212Q_WT < L2FC3,]
mut.rep.both <- target.enhancers[ target.enhancers$padj_E379K_WT < FDR2 & target.enhancers$log2FoldChange_E379K_WT < L2FC3 &
                                    target.enhancers$padj_R212Q_WT < FDR2 & target.enhancers$log2FoldChange_R212Q_WT < L2FC3,]
E379Konly <- mut.rep[ mut.rep$padj_E379K_WT < FDR2 & mut.rep$log2FoldChange_E379K_WT < L2FC3 &
                        !mut.rep$padj_R212Q_WT < FDR2 ,]
R212Qonly <- mut.rep[ !mut.rep$padj_E379K_WT < FDR2 &
                        mut.rep$padj_R212Q_WT < FDR2 & mut.rep$log2FoldChange_R212Q_WT < L2FC3,]

mut.enh <- target.enhancers[ target.enhancers$padj_E379K_WT < FDR2 & target.enhancers$log2FoldChange_E379K_WT > L2FC3 |
                                    target.enhancers$padj_R212Q_WT < FDR2 & target.enhancers$log2FoldChange_R212Q_WT > L2FC3,]
mut.enh.both <- target.enhancers[ target.enhancers$padj_E379K_WT < FDR2 & target.enhancers$log2FoldChange_E379K_WT > L2FC3 &
                               target.enhancers$padj_R212Q_WT < FDR2 & target.enhancers$log2FoldChange_R212Q_WT > L2FC3,]
E379Konly.enh <- mut.enh[ mut.enh$padj_E379K_WT < FDR2 & mut.enh$log2FoldChange_E379K_WT > L2FC3 &
                         !mut.enh$padj_R212Q_WT < FDR2 ,]
R212Qonly.enh <- mut.enh[ mut.enh$padj_R212Q_WT < FDR2 & mut.enh$log2FoldChange_R212Q_WT > L2FC3 &
                            !mut.enh$padj_E379K_WT < FDR2 ,]

non.changing.both <- target.enhancers[target.enhancers$padj_E379K_WT >0.1 & abs(target.enhancers$log2FoldChange_E379K_WT) < 0.25 & target.enhancers$padj_R212Q_WT >0.1 & abs(target.enhancers$log2FoldChange_R212Q_WT) < 0.25, ]
non_sig <- target.enhancers[(!target.enhancers$PeakID %in% mut.rep$PeakID) & !(target.enhancers$PeakID %in% R212Q_induced$PeakID) & !(target.enhancers$PeakID %in% E379K_induced$PeakID),]

#### Fig 7.a ####
a <- ggplot(E379K_constant, aes(x=log(WT_average_norm,10), y=log2FoldChange_E379K_WT))+
        geom_point(size =3, alpha=0.6, shape=16, colour="dark grey")+
        ggtitle("E379K-sensitivity at PPARg2 target enhancers \n1573 enhancers") + 
        ylim(-9,9)+
        geom_point(data=E379K_repressed, aes(x=log(WT_average_norm,10), y=log2FoldChange_E379K_WT), colour="green4", size =3, alpha=0.6, shape=16)+
        geom_point(data=E379K_induced, aes(x=log(WT_average_norm,10), y=log2FoldChange_E379K_WT), colour="red", size =3, alpha=0.6, shape=16)+
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
        geom_hline(yintercept=0)

print(paste0(nrow(E379K_repressed), " peaks gain less PPARg occupancy upon E379K mutation"))
print(paste0(nrow(E379K_induced), " peaks gain more PPARg occupancy upon E379K mutation"))

b <- ggplot(R212Q_constant, aes(x=log(WT_average_norm,10), y=log2FoldChange_R212Q_WT))+
        geom_point(size =3, alpha=0.6, shape=16, colour="dark grey")+
        ggtitle("R212Q-sensitivity at PPARg2 target enhancers \n1573 enhancers") + 
        ylim(-9,9)+
        geom_point(data=R212Q_repressed, aes(x=log(WT_average_norm,10), y=log2FoldChange_R212Q_WT), colour="blue", size =3, alpha=0.6, shape=16)+
        geom_point(data=R212Q_induced, aes(x=log(WT_average_norm,10), y=log2FoldChange_R212Q_WT), colour="red", size =3, alpha=0.6, shape=16)+
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
        geom_hline(yintercept=0)

print(paste0(nrow(R212Q_repressed), " peaks gain less PPARg occupancy upon R212Q mutation"))
print(paste0(nrow(R212Q_induced), " peaks gain more PPARg occupancy upon R212Q mutation"))

pdf('Fig7a_mut_sensitivety.pdf', width=10, height = 5)
grid.arrange(a, b, ncol = 2, nrow = 1)
dev.off()
rm(a,b)

#### Fig 7.b ####
tmp <- target.enhancers[c('chr3-21','chr3-476','chr3-138','chr6-179','chr6-18','chr6-343','chr8-357','chr8-649'),] #selected enhancers

pdf('Fig7b_Mutation_sensitivity_correlation.pdf', width=5, height = 5)
ggplot(non_sig, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT))+
        geom_abline(intercept = 0, slope = 1, linetype="dotted")+
        geom_hline(yintercept = 0)+
        geom_vline(xintercept = 0)+
        geom_point(colour="darkgrey", size=3, shape=16, alpha=0.6) +
        ylim(-8.5,2)+
        xlim(-8.5, 2)+
        geom_point(data=E379Konly, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT), colour="green4", size=3, shape=16, alpha=0.6)+
        geom_point(data=R212Qonly, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT), colour="blue", size=3, shape=16, alpha=0.6)+
        geom_point(data=mut.rep.both, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT), colour="turquoise3", size=3, shape=16, alpha=0.6)+
        geom_point(data=E379Konly.enh, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT), colour="orangered", size=3, shape=16, alpha=0.7)+
        geom_point(data=R212Qonly.enh, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT), colour="orange", size=3, shape=16, alpha=0.7)+
        geom_point(data=mut.enh.both, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT), colour="purple", size=3, shape=16, alpha=0.7)+
        geom_point(data=tmp, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT), colour="black", size=3, shape=1, alpha=1)+
        geom_text_repel(data=tmp, aes(log2FoldChange_E379K_WT, log2FoldChange_R212Q_WT, label=Position))+
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)
dev.off()
rm(tmp)

print(paste0(nrow(E379Konly), " peaks gain less PPARg occupancy upon E379K mutation only"))
print(paste0(nrow(R212Qonly), " peaks gain less PPARg occupancy upon R212Q mutation only"))
print(paste0(nrow(mut.rep.both), " peaks gain less PPARg occupancy upon both E379K and R212Q mutations"))
# Venn diagram plotted using http://eulerr.co/

#### Fig 7.c ####
write.table(non_sig[,c("PeakID", "Chr", "Start", "End", "Strand")], col.names = T, row.names = F, sep="\t", quote=F, file="Insens_enhancers.FACT.txt") 
write.table(mut.rep.both[,c("PeakID", "Chr", "Start", "End", "Strand")], col.names = T, row.names = F, sep="\t", quote=F, file="Dual_enhancers.FACT.txt") 
write.table(E379Konly[,c("PeakID", "Chr", "Start", "End", "Strand")], col.names = T, row.names = F, sep="\t", quote=F, file="E379Konly_enhancers.FACT.txt") 
write.table(R212Qonly[,c("PeakID", "Chr", "Start", "End", "Strand")], col.names = T, row.names = F, sep="\t", quote=F, file="R212Qonly_enhancers.FACT.txt") 

# See separate script for further analysis: BS_enrichment_mult_thresholds.sh


#### Fig 7.e ####
pdf('Fig7e.pdf', width=5, height = 5)
par(pty="s")
boxplot(non_sig$DeNovo_PPRE_TE,
        E379Konly$DeNovo_PPRE_TE,
        R212Qonly$DeNovo_PPRE_TE,
        mut.rep.both$DeNovo_PPRE_TE,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Dual"),
        main = 'DeNovo PPRE (TE) motif score', ylab='Log odds motif score',
        ylim = c(-1,16),
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))
dev.off()

print(paste0(nrow(non_sig), " peaks are insensitive to mutations"))
print(paste0(nrow(E379Konly), " peaks are sensitive to E379K only"))
print(paste0(nrow(R212Qonly), " peaks are sensitive to R212Q only"))
print(paste0(nrow(mut.rep.both), " peaks are sensitive to both mutations"))

#### Fig 7.f ####
pdf('Fig7f.pdf', width=10, height = 5)
par(mfrow=c(1,2), pty="s")
boxplot(non_sig$MED1_Control_average_norm,
        E379Konly$MED1_Control_average_norm,
        R212Qonly$MED1_Control_average_norm,
        mut.rep.both$MED1_Control_average_norm,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Dual"),
        main = 'MED1 ChIP-seq signal, basal level', ylab='MED1 Control_average_norm',
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))

boxplot(non_sig$H3K27ac_Control_average_norm,
        E379Konly$H3K27ac_Control_average_norm,
        R212Qonly$H3K27ac_Control_average_norm,
        mut.rep.both$H3K27ac_Control_average_norm,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Dual"),
        main = 'H3K27ac ChIP-seq signal, basal level', ylab='H3K27ac_Control_average_norm',
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))
dev.off()

#### Fig 7.g ####
pdf('Fig7g.pdf', width=10, height = 5)
par(mfrow=c(1,2), pty="s")
boxplot(non_sig$MED1_log2FoldChange_WT_Control,
        E379Konly$MED1_log2FoldChange_WT_Control,
        R212Qonly$MED1_log2FoldChange_WT_Control,
        mut.rep.both$MED1_log2FoldChange_WT_Control,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Dual"),
        main = 'MED1 ChIP-seq signal, PPARg-WT', ylab='MED1 log2FoldChange_WT_Control',
        ylim = c(-0.25,5.7),
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))

boxplot(non_sig$H3K27ac_log2FoldChange_WT_Control,
        E379Konly$H3K27ac_log2FoldChange_WT_Control,
        R212Qonly$H3K27ac_log2FoldChange_WT_Control,
        mut.rep.both$H3K27ac_log2FoldChange_WT_Control,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Dual"),
        main = 'H3K27ac ChIP-seq signal, PPARg-WT', ylab='H3K27ac log2FoldChange_WT_Control',
        ylim = c(-0.15,3.7),
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))
dev.off()

#### Statistical test Fig 7.e-g ####
# setup dataframe for pairwise statistical analysis:
mut.rep.both$EnhancerGroup <- 'Dual'
non_sig$EnhancerGroup <- 'Insens'
E379Konly$EnhancerGroup <- 'E379Konly'
R212Qonly$EnhancerGroup <- 'R212Qonly'
d.mut.effect <- rbind(non_sig, E379Konly, R212Qonly, mut.rep.both)

# pairwise wilcoxon rank sum test:
pairwise.wilcox.test(d.mut.effect$DeNovo_PPRE_TE, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.mut.effect$MED1_Control_average_norm, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.mut.effect$H3K27ac_Control_average_norm, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.mut.effect$MED1_log2FoldChange_WT_Control, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.mut.effect$H3K27ac_log2FoldChange_WT_Control, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")




# _______________________________________  Figure 8 __________________________________________####

#### Integrating data ####
# read in ATAC-aggregation data and devide in groups:
Agg.ATAC <- read.delim("TagCounts_WTpeaksTH35Ext_ATAC_hist_adj_mm10.txt")
# Order of samples:
# ATAC_Control_pool.TD [,2:62]
# ATAC_E379K_pool.TD [,63:123]
# ATAC_R212Q_pool.TD [,124:184]
# ATAC_WT_pool.TD [,185:245]
# ControlRep1_star.Aligned.out.primary.dedup_under120.TD [,246:306]
# ControlRep2_star.Aligned.out.primary.dedup_under120.TD [,307:367]
# E379KRep1_star.Aligned.out.primary.dedup_under120.TD [,368:428]
# E379KRep2_star.Aligned.out.primary.dedup_under120.TD [,429:489]
# R212QRep1_star.Aligned.out.primary.dedup_under120.TD [,490:550]
# R212QRep2_star.Aligned.out.primary.dedup_under120.TD [,551:611]
# WTRep1_star.Aligned.out.primary.dedup_under120.TD [,612:672]
# WTRep2_star.Aligned.out.primary.dedup_under120.TD [,673:733]
colnames(Agg.ATAC)[1] <- 'PeakID'

ATAC.hist.mut.rep <- mut.rep[,1:2]
ATAC.hist.mut.rep <- merge(ATAC.hist.mut.rep, Agg.ATAC, by='PeakID')
ATAC.hist.mut.rep$Chr <- NULL

ATAC.hist.E379Konly <- E379Konly[,1:2]
ATAC.hist.E379Konly <- merge(ATAC.hist.E379Konly, Agg.ATAC, by='PeakID')
ATAC.hist.E379Konly$Chr <- NULL

ATAC.hist.R212Qonly <- R212Qonly[,1:2]
ATAC.hist.R212Qonly <- merge(ATAC.hist.R212Qonly, Agg.ATAC, by='PeakID')
ATAC.hist.R212Qonly$Chr <- NULL

ATAC.hist.mut.rep.both <- mut.rep.both[,1:2]
ATAC.hist.mut.rep.both <- merge(ATAC.hist.mut.rep.both, Agg.ATAC, by='PeakID')
ATAC.hist.mut.rep.both$Chr <- NULL

ATAC.hist.non_sig <- non_sig[,1:2]
ATAC.hist.non_sig <- merge(ATAC.hist.non_sig, Agg.ATAC, by='PeakID')
ATAC.hist.non_sig$Chr <- NULL

ATAC.target.enhancers <- target.enhancers[,1:2]
ATAC.target.enhancers <- merge(ATAC.target.enhancers, Agg.ATAC, by='PeakID')
ATAC.target.enhancers$Chr <- NULL

# Setting up the data-frame for aggregation-plot using ggplot. The new data-frame contains the mean across all enhancers (in the selected group) at each position relative to peak-center. 
# Target enhancers
dat1 <- as.data.frame(rep("Control",61)) # generates a data-frame where all rows have a column with the id 'Control'. It contain 61 rows corresponding to number of bins in the ATAC-aggregation count matrix.
dat1$rep1 <- apply(ATAC.target.enhancers[,246:306], 2, mean) # this takes the mean across all enhancers in a bin, through all bins
dat1$rep2 <- apply(ATAC.target.enhancers[,307:367], 2, mean) # as above
dat1$comb <- apply(ATAC.target.enhancers[,2:62], 2, mean) # as above
dat1$bin <- seq(1:61) # makes a new column named 'bin' and numbers from 1 to 61
colnames(dat1)[1] <- "id" # renames column 1 to 'id'

dat2 <- as.data.frame(rep("WT",61))
dat2$rep1 <- apply(ATAC.target.enhancers[,612:672], 2, mean)
dat2$rep2 <- apply(ATAC.target.enhancers[,673:733], 2, mean)
dat2$comb <- apply(ATAC.target.enhancers[,185:245], 2, mean)
dat2$bin <- seq(1:61)
colnames(dat2)[1] <- "id"

dat3 <- as.data.frame(rep("E379K",61))
dat3$rep1 <- apply(ATAC.target.enhancers[,368:428], 2, mean)
dat3$rep2 <- apply(ATAC.target.enhancers[,429:489], 2, mean)
dat3$comb <- apply(ATAC.target.enhancers[,63:123], 2, mean)
dat3$bin <- seq(1:61)
colnames(dat3)[1] <- "id"

dat4 <- as.data.frame(rep("R212Q",61))
dat4$rep1 <- apply(ATAC.target.enhancers[,490:550], 2, mean)
dat4$rep2 <- apply(ATAC.target.enhancers[,551:611], 2, mean)
dat4$comb <- apply(ATAC.target.enhancers[,124:184], 2, mean)
dat4$bin <- seq(1:61)
colnames(dat4)[1] <- "id"

d.ATAC.target.enhancers <- rbind(dat2, dat3, dat4, dat1) # row-bind the four data-frames in the order: 'WT', 'E379K', 'R212Q', 'Control'
d.ATAC.target.enhancers$Enhancer <- "TE" # makes a column named 'Enhancer' and add the identifier 'TE' (target enhancer) to all. 

# Dual sensitive enhancers
dat1 <- as.data.frame(rep("Control",61))
dat1$rep1 <- apply(ATAC.hist.mut.rep.both[,246:306], 2, mean)
dat1$rep2 <- apply(ATAC.hist.mut.rep.both[,307:367], 2, mean)
dat1$comb <- apply(ATAC.hist.mut.rep.both[,2:62], 2, mean)
dat1$bin <- seq(1:61)
colnames(dat1)[1] <- "id"

dat2 <- as.data.frame(rep("WT",61))
dat2$rep1 <- apply(ATAC.hist.mut.rep.both[,612:672], 2, mean)
dat2$rep2 <- apply(ATAC.hist.mut.rep.both[,673:733], 2, mean)
dat2$comb <- apply(ATAC.hist.mut.rep.both[,185:245], 2, mean)
dat2$bin <- seq(1:61)
colnames(dat2)[1] <- "id"

dat3 <- as.data.frame(rep("E379K",61))
dat3$rep1 <- apply(ATAC.hist.mut.rep.both[,368:428], 2, mean)
dat3$rep2 <- apply(ATAC.hist.mut.rep.both[,429:489], 2, mean)
dat3$comb <- apply(ATAC.hist.mut.rep.both[,63:123], 2, mean)
dat3$bin <- seq(1:61)
colnames(dat3)[1] <- "id"

dat4 <- as.data.frame(rep("R212Q",61))
dat4$rep1 <- apply(ATAC.hist.mut.rep.both[,490:550], 2, mean)
dat4$rep2 <- apply(ATAC.hist.mut.rep.both[,551:611], 2, mean)
dat4$comb <- apply(ATAC.hist.mut.rep.both[,124:184], 2, mean)
dat4$bin <- seq(1:61)
colnames(dat4)[1] <- "id"

d.ATAC.dual <- rbind(dat2, dat3, dat4, dat1)
d.ATAC.dual$Enhancer <- "Dual"

# E379Konly sensitive enhancers
dat1 <- as.data.frame(rep("Control",61))
dat1$rep1 <- apply(ATAC.hist.E379Konly[,246:306], 2, mean)
dat1$rep2 <- apply(ATAC.hist.E379Konly[,307:367], 2, mean)
dat1$comb <- apply(ATAC.hist.E379Konly[,2:62], 2, mean)
dat1$bin <- seq(1:61)
colnames(dat1)[1] <- "id"

dat2 <- as.data.frame(rep("WT",61))
dat2$rep1 <- apply(ATAC.hist.E379Konly[,612:672], 2, mean)
dat2$rep2 <- apply(ATAC.hist.E379Konly[,673:733], 2, mean)
dat2$comb <- apply(ATAC.hist.E379Konly[,185:245], 2, mean)
dat2$bin <- seq(1:61)
colnames(dat2)[1] <- "id"

dat3 <- as.data.frame(rep("E379K",61))
dat3$rep1 <- apply(ATAC.hist.E379Konly[,368:428], 2, mean)
dat3$rep2 <- apply(ATAC.hist.E379Konly[,429:489], 2, mean)
dat3$comb <- apply(ATAC.hist.E379Konly[,63:123], 2, mean)
dat3$bin <- seq(1:61)
colnames(dat3)[1] <- "id"

dat4 <- as.data.frame(rep("R212Q",61))
dat4$rep1 <- apply(ATAC.hist.E379Konly[,490:550], 2, mean)
dat4$rep2 <- apply(ATAC.hist.E379Konly[,551:611], 2, mean)
dat4$comb <- apply(ATAC.hist.E379Konly[,124:184], 2, mean)
dat4$bin <- seq(1:61)
colnames(dat4)[1] <- "id"

d.ATAC.E379Konly <- rbind(dat2, dat3, dat4, dat1)
d.ATAC.E379Konly$Enhancer <- "E379Konly"

# R212Qonly sensitive enhancers
dat1 <- as.data.frame(rep("Control",61))
dat1$rep1 <- apply(ATAC.hist.R212Qonly[,246:306], 2, mean)
dat1$rep2 <- apply(ATAC.hist.R212Qonly[,307:367], 2, mean)
dat1$comb <- apply(ATAC.hist.R212Qonly[,2:62], 2, mean)
dat1$bin <- seq(1:61)
colnames(dat1)[1] <- "id"

dat2 <- as.data.frame(rep("WT",61))
dat2$rep1 <- apply(ATAC.hist.R212Qonly[,612:672], 2, mean)
dat2$rep2 <- apply(ATAC.hist.R212Qonly[,673:733], 2, mean)
dat2$comb <- apply(ATAC.hist.R212Qonly[,185:245], 2, mean)
dat2$bin <- seq(1:61)
colnames(dat2)[1] <- "id"

dat3 <- as.data.frame(rep("E379K",61))
dat3$rep1 <- apply(ATAC.hist.R212Qonly[,368:428], 2, mean)
dat3$rep2 <- apply(ATAC.hist.R212Qonly[,429:489], 2, mean)
dat3$comb <- apply(ATAC.hist.R212Qonly[,63:123], 2, mean)
dat3$bin <- seq(1:61)
colnames(dat3)[1] <- "id"

dat4 <- as.data.frame(rep("R212Q",61))
dat4$rep1 <- apply(ATAC.hist.R212Qonly[,490:550], 2, mean)
dat4$rep2 <- apply(ATAC.hist.R212Qonly[,551:611], 2, mean)
dat4$comb <- apply(ATAC.hist.R212Qonly[,124:184], 2, mean)
dat4$bin <- seq(1:61)
colnames(dat4)[1] <- "id"

d.ATAC.R212Qonly <- rbind(dat2, dat3, dat4, dat1)
d.ATAC.R212Qonly$Enhancer <- "R212Qonly"

# Insensitive enhancers:
dat1 <- as.data.frame(rep("Control",61))
dat1$rep1 <- apply(ATAC.hist.non_sig[,246:306], 2, mean)
dat1$rep2 <- apply(ATAC.hist.non_sig[,307:367], 2, mean)
dat1$comb <- apply(ATAC.hist.non_sig[,2:62], 2, mean)
dat1$bin <- seq(1:61)
colnames(dat1)[1] <- "id"

dat2 <- as.data.frame(rep("WT",61))
dat2$rep1 <- apply(ATAC.hist.non_sig[,612:672], 2, mean)
dat2$rep2 <- apply(ATAC.hist.non_sig[,673:733], 2, mean)
dat2$comb <- apply(ATAC.hist.non_sig[,185:245], 2, mean)
dat2$bin <- seq(1:61)
colnames(dat2)[1] <- "id"

dat3 <- as.data.frame(rep("E379K",61))
dat3$rep1 <- apply(ATAC.hist.non_sig[,368:428], 2, mean)
dat3$rep2 <- apply(ATAC.hist.non_sig[,429:489], 2, mean)
dat3$comb <- apply(ATAC.hist.non_sig[,63:123], 2, mean)
dat3$bin <- seq(1:61)
colnames(dat3)[1] <- "id"

dat4 <- as.data.frame(rep("R212Q",61))
dat4$rep1 <- apply(ATAC.hist.non_sig[,490:550], 2, mean)
dat4$rep2 <- apply(ATAC.hist.non_sig[,551:611], 2, mean)
dat4$comb <- apply(ATAC.hist.non_sig[,124:184], 2, mean)
dat4$bin <- seq(1:61)
colnames(dat4)[1] <- "id"

d.ATAC.insens <- rbind(dat2, dat3, dat4, dat1)
d.ATAC.insens$Enhancer <- "Insens"

# All HA-peaks:
dat1 <- as.data.frame(rep("Control",61))
dat1$rep1 <- apply(Agg.ATAC[,246:306], 2, mean)
dat1$rep2 <- apply(Agg.ATAC[,307:367], 2, mean)
dat1$comb <- apply(Agg.ATAC[,2:62], 2, mean)
dat1$bin <- seq(1:61)
colnames(dat1)[1] <- "id"

dat2 <- as.data.frame(rep("WT",61))
dat2$rep1 <- apply(Agg.ATAC[,612:672], 2, mean)
dat2$rep2 <- apply(Agg.ATAC[,673:733], 2, mean)
dat2$comb <- apply(Agg.ATAC[,185:245], 2, mean)
dat2$bin <- seq(1:61)
colnames(dat2)[1] <- "id"

dat3 <- as.data.frame(rep("E379K",61))
dat3$rep1 <- apply(Agg.ATAC[,368:428], 2, mean)
dat3$rep2 <- apply(Agg.ATAC[,429:489], 2, mean)
dat3$comb <- apply(Agg.ATAC[,63:123], 2, mean)
dat3$bin <- seq(1:61)
colnames(dat3)[1] <- "id"

dat4 <- as.data.frame(rep("R212Q",61))
dat4$rep1 <- apply(Agg.ATAC[,490:550], 2, mean)
dat4$rep2 <- apply(Agg.ATAC[,551:611], 2, mean)
dat4$comb <- apply(Agg.ATAC[,124:184], 2, mean)
dat4$bin <- seq(1:61)
colnames(dat4)[1] <- "id"

d.ATAC.All <- rbind(dat2, dat3, dat4, dat1)
d.ATAC.All$Enhancer <- "AllHA"

# Generate a master data-frame containing mean values for all the different enhancer groups for all positions relative to peak-center:
d.ATAC.hist.master <- rbind(d.ATAC.target.enhancers, d.ATAC.insens, d.ATAC.E379Konly, d.ATAC.R212Qonly, d.ATAC.dual, d.ATAC.All)
rm(dat1, dat2, dat3, dat4, d.ATAC.target.enhancers, d.ATAC.insens,d.ATAC.E379Konly, d.ATAC.R212Qonly, d.ATAC.dual, d.ATAC.All)


#### Fig 8.a ####
tmp <- d.ATAC.hist.master[d.ATAC.hist.master$id == "Control" & (d.ATAC.hist.master$Enhancer == "Insens"| d.ATAC.hist.master$Enhancer == "E379Konly" | d.ATAC.hist.master$Enhancer == "R212Qonly" | d.ATAC.hist.master$Enhancer == "Dual"), ] # select only id "Control" = basal expression, selecting Enhancer "Dual", "E379Konly", "R212Qonly" and "Insens" to get the groupwise view.  
mycolors <- c("turquoise3","green4", "grey", "blue")
pdf('Fig8A_Basal_ATAC_histggplot.pdf', width=5, height = 5)
ggplot(tmp, aes(x=bin, y=comb, group=Enhancer, color = Enhancer)) +
      geom_ribbon(aes(ymin=rep1,ymax=rep2, fill=Enhancer),alpha=0.7)+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
      scale_color_manual(values = mycolors) + # to exclude the legend include ', guide = FALSE' in the bracket here and below
      scale_fill_manual(values = mycolors)+
      ggtitle("Group-dependent basal chromatin accessibility")+
      xlab("Distance from peak center (bp)")+
      ylab("ATAC-seq\nTag count")
dev.off()


#### Fig 8.b ####
tmp <- d.ATAC.hist.master[d.ATAC.hist.master$Enhancer == "AllHA" & (d.ATAC.hist.master$id == "WT"| d.ATAC.hist.master$id == "Control"), ]
mycolors2 <- c("grey", "red")
a <- ggplot(tmp, aes(x=bin, y=comb, group=id, color = id)) +
      geom_ribbon(aes(ymin=rep1,ymax=rep2, fill=id),alpha=0.7)+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
      scale_color_manual(values = mycolors2) + # to exclude the legend include ', guide = FALSE' in the bracket here and below
      scale_fill_manual(values = mycolors2)+
      ggtitle("WT-induced chromatin accessibility at all-HA sites")+
      xlab("Distance from peak center (bp)")+
      ylab("ATAC-seq\nTag count")

tmp <- d.ATAC.hist.master[d.ATAC.hist.master$Enhancer == "TE" & (d.ATAC.hist.master$id == "WT"| d.ATAC.hist.master$id == "Control"), ]
b <- ggplot(tmp, aes(x=bin, y=comb, group=id, color = id)) +
      geom_ribbon(aes(ymin=rep1,ymax=rep2, fill=id),alpha=0.7)+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
      scale_color_manual(values = mycolors2) + # to exclude the legend include ', guide = FALSE' in the bracket here and below
      scale_fill_manual(values = mycolors2)+
      ggtitle("WT-induced chromatin accessibility at target enhancers")+
      xlab("Distance from peak center (bp)")+
      ylab("ATAC-seq\nTag count")


pdf('Fig8B_ATAC_histo_WT.pdf', width=10, height = 5)
grid.arrange(a,b, ncol = 2, nrow = 1)
dev.off()
rm(a,b)

#### Fig 8.c ####
# Defining enhancers that gain/loose accessibility upon PPARg2-WT expression
ATAC.gain <- target.enhancers[target.enhancers$ATAC_padj_WT_Control < FDR2 & target.enhancers$ATAC_log2FoldChange_WT_Control > L2FC2,]
ATAC.loss <- target.enhancers[target.enhancers$ATAC_padj_WT_Control < FDR2 & target.enhancers$ATAC_log2FoldChange_WT_Control < L2FC2,]
ATAC.const <- target.enhancers[!(target.enhancers$PeakID %in% ATAC.gain$PeakID) & !(target.enhancers$PeakID %in% ATAC.loss$PeakID),]

tmp <- target.enhancers[c('chr3-21','chr3-476','chr3-138','chr6-179','chr6-18','chr6-343','chr8-357','chr8-649'),] # selected enhancers

pdf('Fig8C.pdf', width=5, height = 5)
ggplot(ATAC.const, aes(x=log(ATAC_Control_average_norm,10), y=ATAC_log2FoldChange_WT_Control))+
      geom_point(size =3, alpha=0.6, shape=16, colour="dark grey")+
      ggtitle("PPARg-WT induced chromatin remodeling\nat PPARg target enhancers") + 
      ylim(-6.5,6.5)+
      geom_point(data=ATAC.gain, aes(x=log(ATAC_Control_average_norm,10), y=ATAC_log2FoldChange_WT_Control), colour="red", size =3, alpha=0.6, shape=16)+
      geom_point(data=tmp, aes(x=log(ATAC_Control_average_norm,10), y=ATAC_log2FoldChange_WT_Control), colour="black", size =3, alpha=0.6, shape=1)+
      geom_text_repel(data=tmp, aes(log(ATAC_Control_average_norm,10), ATAC_log2FoldChange_WT_Control, label=Position))+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
      geom_hline(yintercept=0)
dev.off()

#### Fig 8.d ####
pdf('Fig8D_boxplot.pdf', width=5, height = 5)
par(pty="s")
boxplot(non_sig$ATAC_log2FoldChange_WT_Control,
        E379Konly$ATAC_log2FoldChange_WT_Control,
        R212Qonly$ATAC_log2FoldChange_WT_Control,
        mut.rep.both$ATAC_log2FoldChange_WT_Control,
        outline=F,
        notch=TRUE,
        ylim = c(-1,4.5),
        main= "ATAC\nEffect of PPARg expression",
        ylab='L2FC(WT vs. control)',
        names=c("Mut-insens.", "E379Konly", "R212Qonly", "Dual sens"),
        las=2,
        col=c("grey","green4","blue", "turquoise3"))
abline(h=0, lty=2)
dev.off()

pairwise.wilcox.test(d.mut.effect$ATAC_log2FoldChange_WT_Control, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")


#### Fig 8.e ####
pdf('Fig8E_boxplot.pdf', width=5, height = 5)
par(pty="s")
boxplot(target.enhancers$ATAC_Control_average_norm,
        target.enhancers$ATAC_WT_average_norm,
        target.enhancers$ATAC_E379K_average_norm,
        target.enhancers$ATAC_R212Q_average_norm,
        outline=F,
        notch=TRUE,
        ylim = c(-1,85),
        main= "ATAC\ntarget enhancers",
        ylab='Normalized counts',
        names=c( "Control","WT", "E379K", "R212Q"),
        las=2,
        col=c("grey","red", "green4", "blue"))
dev.off()

dat1 <- as.data.frame(rep("Control",1573)) # generates a data-frame where all rows have a column with the id 'Control'. It contain 1573 rows corresponding to number of target enhancers.
dat1$ATAC_normalized_count <- target.enhancers$ATAC_Control_average_norm
colnames(dat1)[1] <- "Treatment"
dat2 <- as.data.frame(rep("WT",1573)) 
dat2$ATAC_normalized_count <- target.enhancers$ATAC_WT_average_norm
colnames(dat2)[1] <- "Treatment"
dat3 <- as.data.frame(rep("E379K",1573)) 
dat3$ATAC_normalized_count <- target.enhancers$ATAC_E379K_average_norm
colnames(dat3)[1] <- "Treatment"
dat4 <- as.data.frame(rep("R212Q",1573))
dat4$ATAC_normalized_count <- target.enhancers$ATAC_R212Q_average_norm
colnames(dat4)[1] <- "Treatment"
d.ATAC.treatment <- rbind(dat1,dat2, dat3, dat4)
pairwise.wilcox.test(d.ATAC.treatment$ATAC_normalized_count, d.ATAC.treatment$Treatment, p.adjust.method = "BH")

#### Fig 8.f ####
# Defining WT-remodeled enhancers that gain/loose accessibility dependent on mutation
ATAC.gain_E379Krep <- ATAC.gain[ATAC.gain$ATAC_padj_E379K_WT < FDR2 & ATAC.gain$ATAC_log2FoldChange_E379K_WT< L2FC2,]
ATAC.gain_R212Qrep <- ATAC.gain[ATAC.gain$ATAC_padj_R212Q_WT < FDR2 & ATAC.gain$ATAC_log2FoldChange_R212Q_WT< L2FC2,]

tmp <- target.enhancers[c('chr3-21','chr3-476'),] 
a <- ggplot(ATAC.gain, aes(x=ATAC_WT_average_norm, y=ATAC_E379K_average_norm))+
      geom_point(size =3, alpha=0.6, shape=16, colour="light green")+
      ggtitle("E379K vs. WT") + 
      geom_abline(intercept = 0, slope = 1, linetype="dotted")+
      ylim(0,200)+
      xlim(0,200)+
      geom_point(data=ATAC.gain_E379Krep, aes(x=ATAC_WT_average_norm, y=ATAC_E379K_average_norm), colour="green4", size =3, alpha=0.6, shape=16)+
      geom_point(data=tmp, aes(x=ATAC_WT_average_norm, y=ATAC_E379K_average_norm), colour="black", size =3, alpha=1, shape=1)+
      geom_text_repel(data=tmp, aes(ATAC_WT_average_norm, ATAC_E379K_average_norm, label=Position))+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
      geom_smooth(data=ATAC.gain, method = "lm", se=FALSE, colour='black', linetype='dashed')

b <- ggplot(ATAC.gain, aes(x=ATAC_WT_average_norm, y=ATAC_R212Q_average_norm))+
      geom_point(size =3, alpha=0.6, shape=16, colour="light blue")+
      ggtitle("R212Q vs. WT") + 
      geom_abline(intercept = 0, slope = 1, linetype="dotted")+
      ylim(0,200)+
      xlim(0,200)+
      geom_point(data=ATAC.gain_R212Qrep, aes(x=ATAC_WT_average_norm, y=ATAC_R212Q_average_norm), colour="blue", size =3, alpha=0.6, shape=16)+
      geom_point(data=tmp, aes(x=ATAC_WT_average_norm, y=ATAC_R212Q_average_norm), colour="black", size =3, alpha=1, shape=1)+
      geom_text_repel(data=tmp, aes(ATAC_WT_average_norm, ATAC_R212Q_average_norm, label=Position))+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)+
      geom_smooth(data=ATAC.gain, method = "lm", se=FALSE, colour='black', linetype='dashed')

pdf('Fig8F_mutation effect.pdf', width=10, height = 5)
grid.arrange(a,b, ncol = 2, nrow = 1)
dev.off()
rm(a,b)

#### Fig 8.g ####
pdf('Fig8G_boxplot.pdf', width=5, height = 5)
par(mfrow=c(1,2), pty="s")
boxplot(non_sig$ATAC_log2FoldChange_E379K_WT,
        E379Konly$ATAC_log2FoldChange_E379K_WT,
        R212Qonly$ATAC_log2FoldChange_E379K_WT,
        mut.rep.both$ATAC_log2FoldChange_E379K_WT,
        outline=F,
        notch=TRUE,
        main= "ATAC\nEffect of E379K mutation",
        ylab='L2FC(Mut vs. WT)',
        ylim=c(-3,1.5),
        names=c("Mut-insens.", "E379Konly", "R212Qonly", "Dual sens"),
        las=2,
        col=c("grey","green4","blue", "turquoise3"))
abline(h=0, lty=2)

boxplot(non_sig$ATAC_log2FoldChange_R212Q_WT,
        E379Konly$ATAC_log2FoldChange_R212Q_WT,
        R212Qonly$ATAC_log2FoldChange_R212Q_WT,
        mut.rep.both$ATAC_log2FoldChange_R212Q_WT,
        outline=F,
        notch=TRUE,
        main= "ATAC\nEffect of R212Q mutation",
        ylab='L2FC(Mut vs. WT)',
        ylim=c(-3,1.5),
        names=c("Mut-insens.", "E379Konly", "R212Qonly", "Dual sens"),
        las=2,
        col=c("grey","green4","blue", "turquoise3"))
abline(h=0, lty=2)
dev.off()

pairwise.wilcox.test(d.mut.effect$ATAC_log2FoldChange_E379K_WT, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.mut.effect$ATAC_log2FoldChange_R212Q_WT, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")

#### Fig 8.h ####
# Integrating data from DOI: 10.1128/MCB.01344-13
d.PPARgCEBPa <- read.delim("HA_peaks_35_PPARgCEBPa_syn.txt")
d.PPARgCEBPa
d.chip.master <- merge(d.chip.master, d.PPARgCEBPa, by="PeakID")
rownames(d.chip.master) <- d.chip.master$PeakID

target.enhancers <- merge(target.enhancers, d.PPARgCEBPa, by="PeakID")
rownames(target.enhancers) <- target.enhancers$PeakID

non_sig <- merge(non_sig, d.PPARgCEBPa, by="PeakID")
rownames(non_sig) <- non_sig$PeakID

mut.rep.both <- merge(mut.rep.both, d.PPARgCEBPa, by="PeakID")
rownames(mut.rep.both) <- mut.rep.both$PeakID

E379Konly <- merge(E379Konly, d.PPARgCEBPa, by="PeakID")
rownames(E379Konly) <- E379Konly$PeakID

R212Qonly <- merge(R212Qonly, d.PPARgCEBPa, by="PeakID")
rownames(R212Qonly) <- R212Qonly$PeakID

# In each enhancer group find enhancers where C/EBPa is bound (CEBPa_bind ==1) and where PPARg facilitates gain in C/EBPa binding (PPARg_eff=='Gain'):
tmp <- non_sig[non_sig$CEBPa_bind==1,] # 277 of 519
tmp1 <- tmp[tmp$PPARg_eff=='Gain',] # 83
tmp <- E379Konly[E379Konly$CEBPa_bind==1,] # 337 of 520
tmp2 <- tmp[tmp$PPARg_eff=='Gain',] # 92
tmp <- R212Qonly[R212Qonly$CEBPa_bind==1,] # 37 of 109
tmp3 <- tmp[tmp$PPARg_eff=='Gain',] # 20
tmp <- mut.rep.both[mut.rep.both$CEBPa_bind==1,] # 202 of 405
tmp4 <- tmp[tmp$PPARg_eff=='Gain',] # 146
tmp <- rbind(tmp1, tmp2, tmp3, tmp4)

# boxplot
pdf('Fig8H_Synergy_PPARgCEBPa_enhancer_groups.pdf', width=5, height = 5)
par(mfrow=c(1,1), pty="s")
boxplot(log2(tmp1$SC14AA_PPARgCEBPa/tmp1$SC14AA_CEBPa),
        log2(tmp2$SC14AA_PPARgCEBPa/tmp2$SC14AA_CEBPa),
        log2(tmp3$SC14AA_PPARgCEBPa/tmp3$SC14AA_CEBPa),
        log2(tmp4$SC14AA_PPARgCEBPa/tmp4$SC14AA_CEBPa),
        outline=F,
        notch=T,
        ylim=c(-0.5, 7.5),
        names=c("insens", "E379K only", "R212Q only", "dual sens"),
        main = 'PPARg-facilitated change in CEBPa binding at gained sites', ylab='C/EBPa ChIP-seq log2(PPARgCEBPa/CEBPa)',
        las=2,
        col=c('grey', "green4", "blue",'turquoise3'))
abline(h=0, lty=2)
dev.off()

pairwise.wilcox.test(log2(tmp$SC14AA_PPARgCEBPa/tmp$SC14AA_CEBPa), tmp$EnhancerGroup, p.adjust.method = "BH")


# _______________________________________  Figure 9 __________________________________________####
#### Dissect motif score ####
# Dissection of the motif strengt within the motif. First part of the analysis is run on the server:
# Extract the sequence of each HA-peak that has passed the threshold, 200 bp at peak center. Run:   ./extractSequence HA_peaks_35_200bp.pos /data/Genomes/mouse/mm10/mm10.fa > sequences_200bp.bed
# Manually change name sequences_200bp.bed > sequences_200bp.tab
# Search for the deNovoMotif within the extracted sequences using the deNovo top scoring motif. The threshold for calling a motif is set down to -2. Run:   homer2 find -s sequences_200bp.tab -m motifRV1_THm2.motif > motif_200bp_THm2.scores

d.Motif <- read.delim("motif_200bp_THm2.scores", header=F)
PWMforward <- read.delim("motif1RV.motif", skip=1, header=F) # we consider the RV motif as the forward
PWMreverse <- read.delim("motif1.motif", skip=1, header=F) # we consider the FW motif as the reverse
colnames(PWMforward) <- c("A","C","G","T")
colnames(PWMreverse) <- c("A","C","G","T")

d.motif <- d.motif[ order(d.motif$V1, -d.motif$V6),]  # this will order motifs with same PeakID with the strongest motif first
d.motif <- d.motif[ duplicated(d.motif$V1)==F,] # this will only allow 1 line per PeakID. -> as the motifs are ordered above, the strongest motif will be kept.

d.motif$Score <- 0
for (i in 1:nrow(d.motif)) {
  Score <- 0
  if (d.motif[i,5] == "+") {
    for (q in 1:17) {
      if (substr(as.character(d.motif[i,3]),q,q) == "A") {
        Score <- Score + log(PWMforward[q,1] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "C") {
        Score <- Score + log(PWMforward[q,2] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "G") {
        Score <- Score + log(PWMforward[q,3] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "T") {
        Score <- Score + log(PWMforward[q,4] / 0.25)
      }
    }
  } else {
    for (q in 1:17) {
      if (substr(as.character(d.motif[i,3]),q,q) == "A") {
        Score <- Score + log(PWMreverse[q,1] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "C") {
        Score <- Score + log(PWMreverse[q,2] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "G") {
        Score <- Score + log(PWMreverse[q,3] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "T") {
        Score <- Score + log(PWMreverse[q,4] / 0.25)
      }
    }
  }
  d.motif[i,"Score"] <- Score
}

d.motif$FivePrime_Score <- 0
for (i in 1:nrow(d.motif)) {
  Score <- 0
  if (d.motif[i,5] == "+") {
    for (q in 1:4) {
      if (substr(as.character(d.motif[i,3]),q,q) == "A") {
        Score <- Score + log(PWMforward[q,1] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "C") {
        Score <- Score + log(PWMforward[q,2] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "G") {
        Score <- Score + log(PWMforward[q,3] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "T") {
        Score <- Score + log(PWMforward[q,4] / 0.25)
      }
    }
  } else {
    for (q in 14:17) {
      if (substr(as.character(d.motif[i,3]),q,q) == "A") {
        Score <- Score + log(PWMreverse[q,1] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "C") {
        Score <- Score + log(PWMreverse[q,2] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "G") {
        Score <- Score + log(PWMreverse[q,3] / 0.25)
      }
      if (substr(as.character(d.motif[i,3]),q,q) == "T") {
        Score <- Score + log(PWMreverse[q,4] / 0.25)
      }
    }
  }
d.motif[i,"FivePrime_Score"] <- Score


d.motif$PPAR_HS_Score <- 0
Score <- 0
if (d.motif[i,5] == "+") {
  for (q in 5:10) {
    if (substr(as.character(d.motif[i,3]),q,q) == "A") {
      Score <- Score + log(PWMforward[q,1] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "C") {
      Score <- Score + log(PWMforward[q,2] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "G") {
      Score <- Score + log(PWMforward[q,3] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "T") {
      Score <- Score + log(PWMforward[q,4] / 0.25)
    }
  }
} else {
  for (q in 8:13) {
    if (substr(as.character(d.motif[i,3]),q,q) == "A") {
      Score <- Score + log(PWMreverse[q,1] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "C") {
      Score <- Score + log(PWMreverse[q,2] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "G") {
      Score <- Score + log(PWMreverse[q,3] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "T") {
      Score <- Score + log(PWMreverse[q,4] / 0.25)
    }
  }
}
d.motif[i,"PPAR_HS_Score"] <- Score


d.motif$RXR_HS_Score <- 0
Score <- 0
if (d.motif[i,5] == "+") {
  for (q in 12:17) {
    if (substr(as.character(d.motif[i,3]),q,q) == "A") {
      Score <- Score + log(PWMforward[q,1] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "C") {
      Score <- Score + log(PWMforward[q,2] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "G") {
      Score <- Score + log(PWMforward[q,3] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "T") {
      Score <- Score + log(PWMforward[q,4] / 0.25)
    }
  }
} else {
  for (q in 1:6) {
    if (substr(as.character(d.motif[i,3]),q,q) == "A") {
      Score <- Score + log(PWMreverse[q,1] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "C") {
      Score <- Score + log(PWMreverse[q,2] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "G") {
      Score <- Score + log(PWMreverse[q,3] / 0.25)
    }
    if (substr(as.character(d.motif[i,3]),q,q) == "T") {
      Score <- Score + log(PWMreverse[q,4] / 0.25)
    }
  }
}
d.motif[i,"RXR_HS_Score"] <- Score
print(i)
}
  
#### Integrate data ####
colnames(d.motif)[1] <- 'PeakID'
# merging groups into motif data frame
d.Motif.mut.rep.both <- mut.rep.both[,1:2]
d.Motif.mut.rep.both <- merge(d.Motif.mut.rep.both, d.motif, by='PeakID')
d.Motif.mut.rep.both$Chr <- NULL
d.Motif.mut.rep.both$EnhancerGroup <- 'Dual'

d.Motif.E379Konly <- E379Konly[,1:2]
d.Motif.E379Konly <- merge(d.Motif.E379Konly, d.motif, by='PeakID')
d.Motif.E379Konly$Chr <- NULL
d.Motif.E379Konly$EnhancerGroup <- 'E379Konly'

d.Motif.R212Qonly <- R212Qonly[,1:2]
d.Motif.R212Qonly <- merge(d.Motif.R212Qonly, d.motif, by='PeakID')
d.Motif.R212Qonly$Chr <- NULL
d.Motif.R212Qonly$EnhancerGroup <- 'R212Qonly'

d.Motif.non_sig <- non_sig[,1:2]
d.Motif.non_sig <- merge(d.Motif.non_sig, d.motif, by='PeakID')
d.Motif.non_sig$Chr <- NULL
d.Motif.non_sig$EnhancerGroup <- 'Insens'

d.Motif.master <- rbind(d.Motif.non_sig, d.Motif.E379Konly,  d.Motif.R212Qonly, d.Motif.mut.rep.both)


#### Fig 9.a ####
pdf('9A.pdf', width=15, height = 5)
par(mfrow=c(1,3), pty="s")
boxplot(d.Motif.non_sig$FivePrime_Score,
        d.Motif.E379Konly$FivePrime_Score,
        d.Motif.R212Qonly$FivePrime_Score,
        d.Motif.mut.rep.both$FivePrime_Score,
        outline=F,
        notch=T,
        ylim=c(-4, 8),
        names=c("insens", "E379K only", "R212Q only", "dual sens"),
        main = 'Motif score 5 prime extension', ylab='log odds motif score',
        las=2,
        col=c('grey', "green4", "blue",'turquoise3'))

boxplot(d.Motif.non_sig$PPAR_HS_Score,
        d.Motif.E379Konly$PPAR_HS_Score,
        d.Motif.R212Qonly$PPAR_HS_Score,
        d.Motif.mut.rep.both$PPAR_HS_Score,
        outline=F,
        notch=T,
        ylim=c(-4, 8),
        names=c("insens", "E379K only", "R212Q only", "dual sens"),
        main = 'Motif score PPAR half site', ylab='log odds motif score',
        las=2,
        col=c('grey', "green4", "blue",'turquoise3'))

boxplot(d.Motif.non_sig$RXR_HS_Score,
        d.Motif.E379Konly$RXR_HS_Score,
        d.Motif.R212Qonly$RXR_HS_Score,
        d.Motif.mut.rep.both$RXR_HS_Score,
        outline=F,
        notch=T,
        ylim=c(-4, 8),
        names=c("insens", "E379K only", "R212Q only", "dual sens"),
        main = 'Motif score RXR half site', ylab='log odds motif score',
        las=2,
        col=c('grey', "green4", "blue",'turquoise3'))
dev.off()

print(paste0(nrow(d.Motif.non_sig), " enhancers are insensitive to mutations"))
print(paste0(nrow(d.Motif.E379Konly), " enhancers are only sensitive to E379K mutations"))
print(paste0(nrow(d.Motif.R212Qonly), " enhancers are only sensitive to R212Q mutations"))
print(paste0(nrow(d.Motif.mut.rep.both), " enhancers are sensitive to both mutations"))

# Pairwise-Wilcox test: 
pairwise.wilcox.test(d.Motif.master$FivePrime_Score, d.Motif.master$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.Motif.master$PPAR_HS_Score, d.Motif.master$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.Motif.master$RXR_HS_Score, d.Motif.master$EnhancerGroup, p.adjust.method = "BH")

#### Defining groups of enhancer accessibility: ####
Open <- 15
Acc.Enhancer <- target.enhancers[target.enhancers$ATAC_Control_average_norm > Open,]
Inacc.Enhancer <- target.enhancers[target.enhancers$ATAC_Control_average_norm < Open,]

Acc.Enhancer_1 <- Acc.Enhancer[Acc.Enhancer$DeNovo_PPRE_TE > 12,]
Acc.Enhancer_2 <- Acc.Enhancer[Acc.Enhancer$DeNovo_PPRE_TE > 11 & Acc.Enhancer$DeNovo_PPRE_TE < 12,]
Acc.Enhancer_3 <- Acc.Enhancer[Acc.Enhancer$DeNovo_PPRE_TE > 10 & Acc.Enhancer$DeNovo_PPRE_TE < 11,]
Acc.Enhancer_4 <- Acc.Enhancer[Acc.Enhancer$DeNovo_PPRE_TE > 9 & Acc.Enhancer$DeNovo_PPRE_TE < 10,]
Acc.Enhancer_5 <- Acc.Enhancer[Acc.Enhancer$DeNovo_PPRE_TE > 8 & Acc.Enhancer$DeNovo_PPRE_TE < 9,]
Acc.Enhancer_6 <- Acc.Enhancer[Acc.Enhancer$DeNovo_PPRE_TE > 7 & Acc.Enhancer$DeNovo_PPRE_TE < 8,]

Inacc.Enhancer_1 <- Inacc.Enhancer[Inacc.Enhancer$DeNovo_PPRE_TE > 12,]
Inacc.Enhancer_2 <- Inacc.Enhancer[Inacc.Enhancer$DeNovo_PPRE_TE > 11 & Inacc.Enhancer$DeNovo_PPRE_TE < 12,]
Inacc.Enhancer_3 <- Inacc.Enhancer[Inacc.Enhancer$DeNovo_PPRE_TE > 10 & Inacc.Enhancer$DeNovo_PPRE_TE < 11,]
Inacc.Enhancer_4 <- Inacc.Enhancer[Inacc.Enhancer$DeNovo_PPRE_TE > 9 & Inacc.Enhancer$DeNovo_PPRE_TE < 10,]
Inacc.Enhancer_5 <- Inacc.Enhancer[Inacc.Enhancer$DeNovo_PPRE_TE > 8 & Inacc.Enhancer$DeNovo_PPRE_TE < 9,]
Inacc.Enhancer_6 <- Inacc.Enhancer[Inacc.Enhancer$DeNovo_PPRE_TE > 7 & Inacc.Enhancer$DeNovo_PPRE_TE < 8,]

#### Fig 9.b ####
pdf('Fig9B.pdf', width=10, height = 5)
par(mfrow=c(1,2), pty="s")
boxplot(Acc.Enhancer_6$log2FoldChange_E379K_WT,
        Acc.Enhancer_5$log2FoldChange_E379K_WT,
        Acc.Enhancer_4$log2FoldChange_E379K_WT,
        Acc.Enhancer_3$log2FoldChange_E379K_WT,
        Acc.Enhancer_2$log2FoldChange_E379K_WT,
        Acc.Enhancer_1$log2FoldChange_E379K_WT,
        Inacc.Enhancer_6$log2FoldChange_E379K_WT,
        Inacc.Enhancer_5$log2FoldChange_E379K_WT,
        Inacc.Enhancer_4$log2FoldChange_E379K_WT,
        Inacc.Enhancer_3$log2FoldChange_E379K_WT,
        Inacc.Enhancer_2$log2FoldChange_E379K_WT,
        Inacc.Enhancer_1$log2FoldChange_E379K_WT,
        outline=F,
        ylim=c(-6, 1.5),
        main = 'E379K-effect PPARg binding', ylab='PPARg ChIP-seq - L2FC(E379K vs. WT)',
        las=2,
        col=c("white","bisque", "coral",'coral2', "red", "darkred", "white","bisque", "coral",'coral2', "red", "darkred"))
abline(h=0, lty=2)

boxplot(Acc.Enhancer_6$log2FoldChange_R212Q_WT,
        Acc.Enhancer_5$log2FoldChange_R212Q_WT,
        Acc.Enhancer_4$log2FoldChange_R212Q_WT,
        Acc.Enhancer_3$log2FoldChange_R212Q_WT,
        Acc.Enhancer_2$log2FoldChange_R212Q_WT,
        Acc.Enhancer_1$log2FoldChange_R212Q_WT,
        Inacc.Enhancer_6$log2FoldChange_R212Q_WT,
        Inacc.Enhancer_5$log2FoldChange_R212Q_WT,
        Inacc.Enhancer_4$log2FoldChange_R212Q_WT,
        Inacc.Enhancer_3$log2FoldChange_R212Q_WT,
        Inacc.Enhancer_2$log2FoldChange_R212Q_WT,
        Inacc.Enhancer_1$log2FoldChange_R212Q_WT,
        outline=F,
        ylim=c(-6, 1.5),
        main = 'R212Q-effect PPARg binding', ylab='PPARg ChIP-seq - L2FC(R212Q vs. WT)',
        las=2,
        col=c("white","bisque", "coral",'coral2', "red", "darkred", "white","bisque", "coral",'coral2', "red", "darkred"))
abline(h=0, lty=2)
dev.off()

print(paste0(nrow(Acc.Enhancer_1), " enhancers are accessible and have a motif score >12"))
print(paste0(nrow(Acc.Enhancer_2), " enhancers are accessible and have a motif score between 11-12"))
print(paste0(nrow(Acc.Enhancer_3), " enhancers are accessible and have a motif score between 10-11"))
print(paste0(nrow(Acc.Enhancer_4), " enhancers are accessible and have a motif score between 9-10"))
print(paste0(nrow(Acc.Enhancer_5), " enhancers are accessible and have a motif score between 8-9"))
print(paste0(nrow(Acc.Enhancer_6), " enhancers are accessible and have a motif score between 7-8"))

print(paste0(nrow(Inacc.Enhancer_1), " enhancers are inaccessible and have a motif score >12"))
print(paste0(nrow(Inacc.Enhancer_2), " enhancers are inaccessible and have a motif score between 11-12"))
print(paste0(nrow(Inacc.Enhancer_3), " enhancers are inaccessible and have a motif score between 10-11"))
print(paste0(nrow(Inacc.Enhancer_4), " enhancers are inaccessible and have a motif score between 9-10"))
print(paste0(nrow(Inacc.Enhancer_5), " enhancers are inaccessible and have a motif score between 8-9"))
print(paste0(nrow(Inacc.Enhancer_6), " enhancers are inaccessible and have a motif score between 7-8"))


#### Fig 9.c ####
pdf('Fig9C.pdf', width=10, height = 5)
par(mfrow=c(1,2), pty="s")
boxplot(Acc.Enhancer_6$ATAC_log2FoldChange_E379K_WT,
        Acc.Enhancer_5$ATAC_log2FoldChange_E379K_WT,
        Acc.Enhancer_4$ATAC_log2FoldChange_E379K_WT,
        Acc.Enhancer_3$ATAC_log2FoldChange_E379K_WT,
        Acc.Enhancer_2$ATAC_log2FoldChange_E379K_WT,
        Acc.Enhancer_1$ATAC_log2FoldChange_E379K_WT,
        Inacc.Enhancer_6$ATAC_log2FoldChange_E379K_WT,
        Inacc.Enhancer_5$ATAC_log2FoldChange_E379K_WT,
        Inacc.Enhancer_4$ATAC_log2FoldChange_E379K_WT,
        Inacc.Enhancer_3$ATAC_log2FoldChange_E379K_WT,
        Inacc.Enhancer_2$ATAC_log2FoldChange_E379K_WT,
        Inacc.Enhancer_1$ATAC_log2FoldChange_E379K_WT,
        outline=F,
        ylim=c(-3.5, 1.5),
        main = 'E379K-effect remodeling', ylab='ATAC-seq - L2FC(E379K vs. WT)',
        las=2,
        col=c("white","bisque", "coral",'coral2', "red", "darkred", "white","bisque", "coral",'coral2', "red", "darkred"))
abline(h=0, lty=2)

boxplot(Acc.Enhancer_6$ATAC_log2FoldChange_R212Q_WT,
        Acc.Enhancer_5$ATAC_log2FoldChange_R212Q_WT,
        Acc.Enhancer_4$ATAC_log2FoldChange_R212Q_WT,
        Acc.Enhancer_3$ATAC_log2FoldChange_R212Q_WT,
        Acc.Enhancer_2$ATAC_log2FoldChange_R212Q_WT,
        Acc.Enhancer_1$ATAC_log2FoldChange_R212Q_WT,
        Inacc.Enhancer_6$ATAC_log2FoldChange_R212Q_WT,
        Inacc.Enhancer_5$ATAC_log2FoldChange_R212Q_WT,
        Inacc.Enhancer_4$ATAC_log2FoldChange_R212Q_WT,
        Inacc.Enhancer_3$ATAC_log2FoldChange_R212Q_WT,
        Inacc.Enhancer_2$ATAC_log2FoldChange_R212Q_WT,
        Inacc.Enhancer_1$ATAC_log2FoldChange_R212Q_WT,
        outline=F,
        ylim=c(-3.5, 1.5),
        main = 'R212Q-effect remodeling', ylab='ATAC-seq - L2FC(R212Q vs. WT)',
        las=2,
        col=c("white","bisque", "coral",'coral2', "red", "darkred", "white","bisque", "coral",'coral2', "red", "darkred"))
abline(h=0, lty=2)
dev.off()

# ____________________________________ Supp. Figure 5 ________________________________________####
d.roc <- target.enhancers
# class1 = E379K-only sensitive
# class2 = R212Q-only sensitive
# class3 = dual-sensitive

d.roc$class1 <- ifelse(d.roc$padj_E379K_WT < 0.05 & d.roc$log2FoldChange_E379K_WT < 0 & 
                         d.roc$padj_R212Q_WT>0.05, 1, 0)
d.roc$class2 <- ifelse(d.roc$padj_E379K_WT > 0.05 & 
                         d.roc$padj_R212Q_WT<0.05 & d.roc$log2FoldChange_R212Q_WT < 0, 1, 0)
d.roc$class3 <- ifelse(d.roc$padj_E379K_WT < 0.05 & d.roc$log2FoldChange_E379K_WT < 0 & 
                         d.roc$padj_R212Q_WT<0.05 & d.roc$log2FoldChange_R212Q_WT < 0, 1, 0)

table(d.roc$class1)
table(d.roc$class2)
table(d.roc$class3)

# ROC curve
pdf('SuppFig5_ROC.pdf', width=15, height = 5)
par(mfrow=c(1,3), pty="s")
roc(d.roc$class1, d.roc$WT_average_norm, plot=T, col="red", main="E379K-only sensitive")
roc(d.roc$class1, d.roc$DeNovo_PPRE_TE, plot=T, add=T, col="darkred")
roc(d.roc$class1, d.roc$ATAC_Control_average_norm, plot=T, add=T, col="black")
roc(d.roc$class1, d.roc$H3K27ac_Control_average_norm, plot=T, add=T, col="grey")
roc(d.roc$class1, d.roc$MED1_Control_average_norm, plot=T, add=T, col="grey34")

roc(d.roc$class2, d.roc$WT_average_norm, plot=T, col="red", main="R212Q-only sensitive")
roc(d.roc$class2, d.roc$DeNovo_PPRE_TE, plot=T, add=T, col="darkred")
roc(d.roc$class2, d.roc$ATAC_Control_average_norm, plot=T, add=T, col="black")
roc(d.roc$class2, d.roc$H3K27ac_Control_average_norm, plot=T, add=T, col="grey")
roc(d.roc$class2, d.roc$MED1_Control_average_norm, plot=T, add=T, col="grey34")

roc(d.roc$class3, d.roc$WT_average_norm, plot=T, col="red", main="Dual sensitive")
roc(d.roc$class3, d.roc$DeNovo_PPRE_TE, plot=T, add=T, col="darkred")
roc(d.roc$class3, d.roc$ATAC_Control_average_norm, plot=T, add=T, col="black")
roc(d.roc$class3, d.roc$H3K27ac_Control_average_norm, plot=T, add=T, col="grey")
roc(d.roc$class3, d.roc$MED1_Control_average_norm, plot=T, add=T, col="grey34")
dev.off()

# AUC
hold <- as.data.frame(matrix(nrow=5,ncol=2))
colnames(hold) <- c("Name", "AUC")
hold[2,2] <- roc(d.roc$class1, d.roc$DeNovo_PPRE_TE)$auc[1]
hold[1,2] <- roc(d.roc$class1, d.roc$WT_average_norm)$auc[1]
hold[5,2] <- roc(d.roc$class1, d.roc$ATAC_Control_average_norm)$auc[1]
hold[3,2] <- roc(d.roc$class1, d.roc$H3K27ac_Control_average_norm)$auc[1]
hold[4,2] <- roc(d.roc$class1, d.roc$MED1_Control_average_norm)$auc[1]
hold[,1] <- c("WT_HA-ChIP_Signal", "DeNovo_MotifScore", "Basal_H3K27ac", "Basal_MED1", "Basal_ATAC")
hold$Name <- factor(hold$Name, levels=hold$Name)

a <- ggplot(hold, aes(Name, AUC))+
  geom_bar(stat="identity")+
  coord_cartesian(ylim=c(0.5,0.75))+
  geom_hline(yintercept = 0.5, linetype="dotted")+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 2)

roc.test(roc(d.roc$class1, d.roc$DeNovo_PPRE_TE),
         roc(d.roc$class1, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class1, d.roc$H3K27ac_Control_average_norm),
         roc(d.roc$class1, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class1, d.roc$MED1_Control_average_norm),
         roc(d.roc$class1, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class1, d.roc$ATAC_Control_average_norm),
         roc(d.roc$class1, d.roc$WT_average_norm))$p.value

hold <- as.data.frame(matrix(nrow=5,ncol=2))
colnames(hold) <- c("Name", "AUC")
hold[2,2] <- roc(d.roc$class2, d.roc$DeNovo_PPRE_TE)$auc[1]
hold[1,2] <- roc(d.roc$class2, d.roc$WT_average_norm)$auc[1]
hold[5,2] <- roc(d.roc$class2, d.roc$ATAC_Control_average_norm)$auc[1]
hold[3,2] <- roc(d.roc$class2, d.roc$H3K27ac_Control_average_norm)$auc[1]
hold[4,2] <- roc(d.roc$class2, d.roc$MED1_Control_average_norm)$auc[1]
hold[,1] <- c("WT_HA-ChIP_Signal", "DeNovo_MotifScore", "Basal_H3K27ac", "Basal_MED1", "Basal_ATAC")
hold$Name <- factor(hold$Name, levels=hold$Name)

b <- ggplot(hold, aes(Name, AUC))+
  geom_bar(stat="identity")+
  coord_cartesian(ylim=c(0.5,0.75))+
  geom_hline(yintercept = 0.5, linetype="dotted")+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 2)

roc.test(roc(d.roc$class2, d.roc$DeNovo_PPRE_TE),
         roc(d.roc$class2, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class2, d.roc$H3K27ac_Control_average_norm),
         roc(d.roc$class2, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class2, d.roc$MED1_Control_average_norm),
         roc(d.roc$class2, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class2, d.roc$ATAC_Control_average_norm),
         roc(d.roc$class2, d.roc$WT_average_norm))$p.value

hold <- as.data.frame(matrix(nrow=5,ncol=2))
colnames(hold) <- c("Name", "AUC")
hold[2,2] <- roc(d.roc$class3, d.roc$DeNovo_PPRE_TE)$auc[1]
hold[1,2] <- roc(d.roc$class3, d.roc$WT_average_norm)$auc[1]
hold[5,2] <- roc(d.roc$class3, d.roc$ATAC_Control_average_norm)$auc[1]
hold[3,2] <- roc(d.roc$class3, d.roc$H3K27ac_Control_average_norm)$auc[1]
hold[4,2] <- roc(d.roc$class3, d.roc$MED1_Control_average_norm)$auc[1]
hold[,1] <- c("WT_HA-ChIP_Signal", "DeNovo_MotifScore", "Basal_H3K27ac", "Basal_MED1", "Basal_ATAC")
hold$Name <- factor(hold$Name, levels=hold$Name)

c <- ggplot(hold, aes(Name, AUC))+
  geom_bar(stat="identity")+
  coord_cartesian(ylim=c(0.5,0.75))+
  geom_hline(yintercept = 0.5, linetype="dotted")+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 2)

roc.test(roc(d.roc$class3, d.roc$DeNovo_PPRE_TE),
         roc(d.roc$class3, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class3, d.roc$H3K27ac_Control_average_norm),
         roc(d.roc$class3, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class3, d.roc$MED1_Control_average_norm),
         roc(d.roc$class3, d.roc$WT_average_norm))$p.value
roc.test(roc(d.roc$class3, d.roc$ATAC_Control_average_norm),
         roc(d.roc$class3, d.roc$WT_average_norm))$p.value

pdf('SuppFig5_AUC.pdf', width=15, height = 5)
grid.arrange(a,b,c, ncol=3, nrow=1)
dev.off()

rm(a,b,c)
rm(hold)

# ____________________________________ Supp. Figure 7 ________________________________________####
pdf('SuppFig7_mutation_effect_MED1_H3K27ac.pdf', width=10, height = 10)
par(mfrow=c(2,2), pty="s")
boxplot(non_sig$H3K27ac_log2FoldChange_E379K_WT,
        E379Konly$H3K27ac_log2FoldChange_E379K_WT,
        R212Qonly$H3K27ac_log2FoldChange_E379K_WT,
        mut.rep.both$H3K27ac_log2FoldChange_E379K_WT,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Both sens."),
        main = 'H3K27ac ChIP-seq signal', ylab='H3K27ac log2FoldChange_E379K_WT',
        ylim = c(-2.5,0.8),
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))
abline(h=0, lty=2)

boxplot(non_sig$H3K27ac_log2FoldChange_R212Q_WT,
        E379Konly$H3K27ac_log2FoldChange_R212Q_WT,
        R212Qonly$H3K27ac_log2FoldChange_R212Q_WT,
        mut.rep.both$H3K27ac_log2FoldChange_R212Q_WT,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Both sens."),
        main = 'H3K27ac ChIP-seq signal', ylab='H3K27ac log2FoldChange_R212Q_WT',
        ylim = c(-2.5,0.8),
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))
abline(h=0, lty=2)

boxplot(non_sig$MED1_log2FoldChange_E379K_WT,
        E379Konly$MED1_log2FoldChange_E379K_WT,
        R212Qonly$MED1_log2FoldChange_E379K_WT,
        mut.rep.both$MED1_log2FoldChange_E379K_WT,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Both sens."),
        main = 'MED1 ChIP-seq signal', ylab='MED1 log2FoldChange_E379K_WT',
        ylim = c(-4,1.5),
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))
abline(h=0, lty=2)

boxplot(non_sig$MED1_log2FoldChange_R212Q_WT,
        E379Konly$MED1_log2FoldChange_R212Q_WT,
        R212Qonly$MED1_log2FoldChange_R212Q_WT,
        mut.rep.both$MED1_log2FoldChange_R212Q_WT,
        outline=F,
        notch=TRUE,
        names=c("Non sens.", "E379K only", "R212Q only", "Both sens."),
        main = 'MED1 ChIP-seq signal', ylab='MED1 log2FoldChange_R212Q_WT',
        ylim = c(-4,1.5),
        las=2,
        col=c('dark grey', 'green4', 'blue', 'turquoise3'))
abline(h=0, lty=2)
dev.off()

pairwise.wilcox.test(d.mut.effect$H3K27ac_log2FoldChange_E379K_WT, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.mut.effect$H3K27ac_log2FoldChange_R212Q_WT, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.mut.effect$MED1_log2FoldChange_E379K_WT, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")
pairwise.wilcox.test(d.mut.effect$MED1_log2FoldChange_R212Q_WT, d.mut.effect$EnhancerGroup, p.adjust.method = "BH")

