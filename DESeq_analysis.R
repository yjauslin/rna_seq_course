library("DESeq2")
library("ggplot2")
library("clusterProfiler")
library("org.Mm.eg.db")
library("enrichplot")

#adjust to preferred project directory
project_path <- "C:/Users/yanni/OneDrive - Universitaet Bern/Master"

#set working directory to project directory
setwd(project_path)

#read in featureCounts data
counts <- read.table("counts_filtered.txt", header = T)

#save readme data to a vector
readme <- c("SRR7821921", "Lung_WT_Case",
            "SRR7821922",	"Lung_WT_Case",
            "SRR7821918",	"Lung_WT_Case",
            "SRR7821919",	"Lung_WT_Case",
            "SRR7821920",	"Lung_WT_Case",
            "SRR7821937",	"Lung_WT_Control",
            "SRR7821938",	"Lung_WT_Control",
            "SRR7821939",	"Lung_WT_Control",
            "SRR7821949",	"Blood_WT_Case",
            "SRR7821950",	"Blood_WT_Case",
            "SRR7821951",	"Blood_WT_Case",
            "SRR7821952",	"Blood_WT_Case",
            "SRR7821953",	"Blood_WT_Case",
            "SRR7821968",	"Blood_WT_Control",
            "SRR7821969",	"Blood_WT_Control",
            "SRR7821970",	"Blood_WT_Control")

# Reshape the vector into a data frame
readme_df <- data.frame(
  SampleID = readme[seq(1, length(readme), by = 2)],
  Description = readme[seq(2, length(readme), by = 2)],
  stringsAsFactors = FALSE
)

# Split the Description column into components and create col_data object
col_data <- transform(readme_df,
                       Tissue = sapply(strsplit(Description, "_"), `[`, 1),
                       Genotype = sapply(strsplit(Description, "_"), `[`, 2),
                       Condition = sapply(strsplit(Description, "_"), `[`, 3))

#get geneID from counts data
geneID <- counts$Geneid
#create matrix of counts without the geneID column
counts <- as.matrix(counts[,-1])
#set geneIDs as row names
rownames(counts) <- geneID
#change rownames to name of the samples
colnames(counts) <- readme_df$SampleID

#create DESEQDataSet Object and execute analysis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~Tissue+Condition)
#fuse tissue and condition together
dds$group <- factor(paste0(dds$Tissue, dds$Condition))
#set group as design of the DESeq-analysis
design(dds) <- ~ group
#analyze dds for differential gene expression
dds <- DESeq(dds)

#remove dependence of the variance on the mean
vst <- vst(dds, blind = T)

#create PCA
pca_data <- plotPCA(vst, intgroup = c('Tissue', 'Condition'), returnData = T)

#save variance for each component in a variable
percentVar <- round(100 * attr(pca_data, "percentVar")) 

#plot pca
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Tissue))+
  geom_point(size = 3)+
  scale_color_brewer(palette = "Set1") +
  theme_bw()+
  theme(axis.title.x = element_text(size = 16),
           axis.title.y = element_text(size = 16),
           axis.text.x = element_text(size = 14),
           axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16) )+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))  

#store results for blood comparison and look at them in different ways
res_blood <- results(dds, contrast=c("group", "BloodControl", "BloodCase"))
head(res_blood)
summary(res_blood)

#store results for lung comparison and look at them in different ways
res_lung <- results(dds, contrast=c("group", "LungControl", "LungCase"))
head(res_lung)
summary(res_lung)


#define cut-off for the p-value
p_cutoff = 0.01

#define cut-off for log2foldchange
logfold_cutoff = 1

#count number of significant genes for blood
sig_genes_blood <- sum(res_blood$padj < p_cutoff, na.rm = T)
print(sig_genes_blood)

#count number of upregulated genes for blood
upregulated_blood <- sum(res_blood$log2FoldChange > logfold_cutoff & res_blood$padj < p_cutoff, na.rm = T)
print(upregulated_blood)

#count number of downregulated genes for blood
downregulated_blood <- sum(res_blood$log2FoldChange < -logfold_cutoff & res_blood$padj < p_cutoff, na.rm = T)
print(downregulated_blood)



#count number of significant genes for lung
sig_genes_lung <- sum(res_lung$padj < p_cutoff, na.rm = T)
print(sig_genes_lung)

#count number of upregulated genes for lung
upregulated_lung <- sum(res_lung$log2FoldChange > logfold_cutoff & res_lung$padj < p_cutoff, na.rm = T)
print(upregulated_lung)

#count number of downregulated genes for lung
downregulated_lung <- sum(res_lung$log2FoldChange < -logfold_cutoff & res_lung$padj < p_cutoff, na.rm = T)
print(downregulated_lung)


#adjust plot visualization to fit two plots next to each other
par(mfrow = c(1,2), mar = c(4,4,3,2), cex = 0.9, cex.axis = 1.1, cex.lab = 1.1, cex.main = 1.2, cex.lab = 1.1)

# Make a basic volcano plot representing every gene
with(res_blood, plot(log2FoldChange, -log10(padj), pch=20, main="Blood: Control vs. Case", xlim=c(-10,15), ylim = c(0, 450)))

# Add colored points for significant genes in blood comparison (padj < p_cutoff(=0.01): blue if logfoldchange<logfold_cutoff(=2), red if logfoldchange>logfold_cutoff(=2))
with(subset(res_blood, padj<p_cutoff), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res_blood, padj<p_cutoff & log2FoldChange < -logfold_cutoff), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res_blood, padj<p_cutoff & log2FoldChange > logfold_cutoff), points(log2FoldChange, -log10(padj), pch=20, col="red"))

#calculate total genes on the plot
total_genes <- length(rownames(res_blood))
# Add the legend
legend("topright", legend=c(paste0("total (", total_genes, " genes)"), paste0("downregulated (", downregulated_blood, " genes)"), paste0("upregulated (", upregulated_blood, " genes)"), paste0("padj < 0.01 (", sig_genes_blood, " genes)")), col=c("black", "blue", "red","orange"), pch=20, cex = 1.1)


#do the same thing for comparison of lung

# Make a basic volcano plot representing every gene
with(res_lung, plot(log2FoldChange, -log10(padj), pch=20, main="Lung: Control vs. Case", xlim=c(-10,15), ylim = c(0, 450)))

# Add colored points for significant genes in lung comparison (padj < p_cutoff(=0.01): blue if logfoldchange<logfold_cutoff(=2), red if logfoldchange>logfold_cutoff(=2))
with(subset(res_lung, padj<p_cutoff), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res_lung, padj<p_cutoff & log2FoldChange < -logfold_cutoff), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res_lung, padj<p_cutoff & log2FoldChange > logfold_cutoff), points(log2FoldChange, -log10(padj), pch=20, col="red"))

# Add the legend
legend("topright", legend=c(paste0("total (", total_genes, " genes)"), paste0("downregulated (", downregulated_lung, " genes)"), paste0("upregulated (", upregulated_lung, " genes)"), paste0("padj < 0.01 (", sig_genes_lung, " genes)")), col=c("black", "blue", "red","orange"), pch=20, cex = 1.1)


#adjust plot visualization to show multiple plots: 2 rows, 2 columns, adjust size of elements of the plot as well as margins for the plot (mar) and axis elements/labels (mgp)
par(mfrow = c(2,2), mar = c(4,3,3,2), cex.axis = 0.65, cex.lab = 0.7, cex.main = 0.7, cex = 1.5, mgp = c(2, 0.5, 0))

#plot differential expression of some of the genes mentioned in initial publication

#lfit1
padj_lfit1_blood <- res_blood['ENSMUSG00000034459', "padj"]
padj_lfit1_lung <- res_lung['ENSMUSG00000034459', "padj"]
plotCounts(dds, gene='ENSMUSG00000034459', intgroup="group", main = paste("Gene: ", "Ifit1", "\nBloodControl vs. Case: padj = ", format(padj_lfit1_blood, scientific = T, digits = 3) ,
                                                                          "\nLungControl vs. Case: padj = ", format(padj_lfit1_lung, scientific = T, digits = 3)))

#oas1a
padj_oas1a_blood <- res_blood['ENSMUSG00000052776', "padj"]
padj_oas1a_lung <- res_lung['ENSMUSG00000052776', "padj"]
plotCounts(dds, gene='ENSMUSG00000052776', intgroup="group", main = paste("Gene: ", "Oas1a", "\nBloodControl vs. Case: padj = ", format(padj_oas1a_blood, scientific = T, digits = 3),
                                                                          "\n LungControl vs. Case: padj = ", format(padj_oas1a_lung, scientific = T, digits = 3)))

#mx1
padj_mx1_blood <- res_blood['ENSMUSG00000000386', "padj"]
padj_mx1_lung <- res_lung['ENSMUSG00000000386', "padj"]
plotCounts(dds, gene='ENSMUSG00000000386', intgroup="group", main = paste("Gene: ", "Mx1", "\nBloodControl vs. Case: padj = ", format(padj_mx1_blood, scientific = T, digits = 3),
                                                                          "\n LungControl vs. Case: padj = ", format(padj_mx1_lung, scientific = T, digits = 3)))

#irf7
padj_Irf7_blood <- res_blood['ENSMUSG00000025498', "padj"]
padj_Irf7_lung <- res_lung['ENSMUSG00000025498', "padj"]
plotCounts(dds, gene='ENSMUSG00000025498', intgroup="group", main = paste("Gene: ", "Irf7,", "\nBloodControl vs. Case: padj = ", format(padj_Irf7_blood, scientific = T, digits = 3),
                                                                          "\n LungControl vs. Case: padj = ", format(padj_Irf7_lung, scientific = T, digits = 3)))


#filter genes which have padj < p_cutoff and absolute logfoldchange > logfold_cutoff
genes_blood <- rownames(res_blood[(res_blood$padj < p_cutoff & !is.na(res_blood$padj)) & (abs(res_blood$log2FoldChange) > logfold_cutoff),]) 

#do overrepresentation analysis for blood
ego_blood <- enrichGO(gene    = genes_blood, #list of genes to be analyzed
                universe      = rownames(res_blood), #background gene set, entire set of genes to consider
                OrgDb         = org.Mm.eg.db, #organism database, should be mouse
                ont           = "BP", #ontology we're interested in, either MF = molecular function, BP = biological process or CC = cellular
                pAdjustMethod = "BH", #method for adjusting p-value, BH = benjamini-hochberg
                pvalueCutoff  = p_cutoff, #sets threshold for unadjusted p-value
                qvalueCutoff  = p_cutoff, #sets threshold for adjusted p-value
                keyType = "ENSEMBL", #defines format for key identifier
                readable      = TRUE) #makes format readable for us

genes_lung <- rownames(res_lung[(res_lung$padj < p_cutoff & !is.na(res_lung$padj)) & (abs(res_lung$log2FoldChange) > logfold_cutoff),]) 

#do overrepresentation analysis for lung
ego_lung <- enrichGO(gene           = genes_lung, #list of genes to be analyzed
                     universe      = rownames(res_lung), #background gene set, entire set of genes to consider
                     OrgDb         = org.Mm.eg.db, #organism database, should be mouse
                     ont           = "BP", #ontology we're interested in, either MF = molecular function, BP = biological process or CC = cellular
                     pAdjustMethod = "BH", #method for adjusting p-value, BH = benjamini-hochberg
                     pvalueCutoff  = p_cutoff, #sets threshold for unadjusted p-value
                     qvalueCutoff  = p_cutoff, #sets threshold for adjusted p-value
                     keyType = "ENSEMBL", #defines format for key identifier
                     readable      = TRUE) #makes format readable for us


#create barplots for overrepresentation analysis
barplot(ego_blood, showCategory=20) + ggtitle("Blood: Control vs. Case") + theme(text = element_text(size = 16))
barplot(ego_lung, showCategory=20) + ggtitle("Lung: Control vs. Case") + theme(text = element_text(size = 16))


                 