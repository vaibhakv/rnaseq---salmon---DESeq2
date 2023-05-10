library(tximport)
library(DESeq2)
library(pheatmap)
library(apeglm)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(biomaRt)
library(plotly)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)


#creating the metadata from SraRunTable which was downloaded from SRA run selector
#modified a little bit afterwards
metadata <- read_csv("data/metadata.csv")%>% 
  select(`Sample Name`, Genotype) %>% 
  arrange(`Sample Name`) 

metadata2 <- read.csv("data/metadata.csv", row.names = 1, header=TRUE) %>% 
  select(Genotype)
  
#--------------------------------------------------

#creating a sample_files character vector which is required for tx2gene option in tximport
#This basically feeds in the count data (quant.sf files) from their respective folders
sample_files <- paste0('data/', pull(metadata, `Sample Name`), '/quant.sf')
names(sample_files) <- pull(metadata, `Sample Name`)

#creating gene map which is required for tx2gene option in tximport 
#This is done in order to import gene level counts instead of the transcript level counts which were obtained by Salmon
gene_map <- read_csv("data/gene_map.csv", col_names = c("enstid", "ensgid"))

#importing counts data from salmon output into the count_data object
count_data <- tximport(files = sample_files,
                       type = 'salmon',
                       tx2gene = gene_map,
                       ignoreTxVersion = TRUE)

#making sure that the column names in count_data = rownames in metadata
all(colnames(count_data$counts) %in% rownames(metadata2))
all(colnames(count_data$counts) == rownames(metadata2))

#creating DEseq dataset
dds <-  DESeqDataSetFromTximport(txi = count_data,
                                 colData = metadata2,
                                 design = ~Genotype)
counts(dds)
#DON'T RUN THIS THE FIRST TIME------second round-------------------------------------

dds1 <-  DESeqDataSetFromTximport(txi = count_data,
                                 colData = metadata2,
                                 design = ~Genotype)
#DON'T RUN THIS THE FIRST TIME------third round-------------------------------------

dds2 <-  DESeqDataSetFromTximport(txi = count_data,
                                  colData = metadata2,
                                  design = ~Genotype)
#---------------------------------------------------------------------------------------------------------------------------
#This step removes all the rows with sum of the counts less than 10
#I don't know if this is the correct way to to go about this. Might change in future
#If you want to reverse this step, just run the DESeqDataSetFromTximport() again and proceed from there and skip this step

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
#---------------------------------second-round----------------------------------
keep <- rowSums(counts(dds1)) >= 10
dds1 <- dds1[keep, ]
#---------------------------------third-round----------------------------------------------------------------------------
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep, ]

#Estimating normalization factors
sizeFactors <- estimateSizeFactors(dds)
normalizationFactors(sizeFactors)

#Comparing between raw counts and normalized counts
counts(dds) %>% head(5)
counts(sizeFactors, normalized= TRUE) %>% head(5)

#rough plotting -------------------------------
boxplot(counts(sizeFactors, normalized= TRUE))

#log transforming the data - 2 types (vst and rlog)
vst = varianceStabilizingTransformation(dds) #
rld <- rlog(dds, blind=T) 

#basic PCA (PC1, PC2)------------------------- 
plotPCA(object = rld, intgroup='Genotype')   #
plotPCA(object = vst, intgroup='Genotype')   #
#----------------------------------------------

#Advanced PCA----------------------------------------------------------------#
# Input is a matrix of log transformed values                                #
pca <- vst %>%                                                               #
  assay() %>%                                                                #
  t() %>%                                                                    #
  prcomp()                                                                   #
                                                                             #
# Create data frame with metadata and PC3 and PC4 values for input to ggplot #
pca_df <- cbind(metadata2, pca$x)                                             #
ggplot(pca_df) +                                                             #
  geom_point(aes(x=PC3, y=PC4, color = Genotype))                            #
#--------------------------------------------------------------------------- #


#basic dendogram-------------------------------
vst %>% 
  assay() %>% 
  t() %>% 
  dist() %>% 
  hclust() %>% 
  plot()
#-----------------------------------------------

#Heatmap (with rlog)----------------------------------------

#-----------------------------------------------------------
rld %>% 
  assay() %>% 
  cor() %>% 
  pheatmap(annotation = metadata2)

#Heatmap (with vst)-----------------------------------------
vst %>% 
  assay() %>% 
  cor() %>% 
  pheatmap(annotation = metadata2)

#running DESeq!!!-------------------------------------------

#releveling the DESeq dataset, setting the reference to KO as we want to compare the other two conditions to this one
# dds$Genotype <- relevel(dds$Genotype, ref = "KO")
#second-round-------------------------
dds1$Genotype <- relevel(dds1$Genotype, ref = "WT")
#third-round-------------------------
dds2$Genotype <- relevel(dds2$Genotype, ref = "Rescue")

#running DESeq -- reusing the 'dds' object which we used for DESeq dataset. Remember this.
# dds <- DESeq(dds)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

#check the result names. See that the reference is set correct or not.
# resultsNames(dds)
resultsNames(dds1)
resultsNames(dds2)

#get the result table for both comparisons providing the name which we got above.
# res1 <- results(dds, name = "Genotype_WT_vs_KO")
# res2 <- results(dds, name = "Genotype_Rescue_vs_KO")
#second-round-------------
res3 <- results(dds1, name = "Genotype_KO_vs_WT")
#third-round-------------
res4 <- results(dds2, name = "Genotype_KO_vs_Rescue")
#check the summary of both result tables
# summary(res1)
# summary(res2)
summary(res3)
summary(res4)

#Checking the fit of the dispersion estimates
# plotDispEsts(dds)
#second-round-----------------------
plotDispEsts(dds1)
#third-round-----------------------
plotDispEsts(dds2)

#--------------------------RES1-------------------------#
#MA plot for the unshrunk log2 fold changes
# plotMA(res1)
#second-round--------------------
plotMA(res3)
#third-round--------------------
plotMA(res4)
#Shrinking the log 2 fold changes
# res1LFC <- lfcShrink(dds, coef = "Genotype_WT_vs_KO", type = "apeglm")
#second-round--------------------
res3LFC <- lfcShrink(dds1, coef = "Genotype_KO_vs_WT", type = "apeglm")
summary(res3LFC)
plotMA(res3LFC)
#third-round--------------------
res4LFC <- lfcShrink(dds2, coef = "Genotype_KO_vs_Rescue", type = "apeglm")
summary(res4LFC)
plotMA(res4LFC)

# #Note that the summary doesn't change after this step, just the log2 fold change values are changed
# summary(res1LFC)
# #MA plot for the shrunken log2 fold changes
# plotMA(res1LFC)
# #-------------------------------------------------------#
# 
# #--------------------------RES2-------------------------#
# #MA plot for the unshrunk log2 fold changes
# plotMA(res2)
# #Shrinking the log 2 fold changes
# res2LFC <- lfcShrink(dds, coef = "Genotype_Rescue_vs_KO", type = "apeglm")
# #Note that the summary doesn't change after this step, just the log2 fold change values are changed
# summary(res2LFC)
# #MA plot for the shrunken log2 fold changes
# plotMA(res2LFC)
# #-------------------------------------------------------#

#------------------------
padj.cutoff <- 0.05

# res1_tb <- res1 %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble()
# 
# sig1 <- res1_tb %>%
#   filter(padj < padj.cutoff)
# 
# res2_tb <- res2 %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble()
# 
# sig2 <- res2_tb %>%
#   filter(padj < padj.cutoff)

res3_tb <- res3 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#this one I changed the thresholds to visualize the heatmap more comprehensively
sig3 <- res3_tb %>%
  filter(padj < padj.cutoff) #& abs(log2FoldChange) > 2.7)

#just the 
sig3_up <- res3_tb %>%
  filter(padj < padj.cutoff & log2FoldChange > 0)

sig3_down <-  res3_tb %>%
  filter(padj < padj.cutoff & log2FoldChange < 0)


res4_tb <- res4 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig4 <- res4_tb %>%
  filter(padj <  padj.cutoff)

sig4_up <- res4_tb %>%
  filter(padj < padj.cutoff & log2FoldChange > 0)

sig4_down <- res4_tb %>%
  filter(padj < padj.cutoff & log2FoldChange < 0)

metadata2_tb <- metadata2 %>% 
                rownames_to_column(var = "samplename") %>% 
                as.tibble()

normalized_counts <- counts(sizeFactors, normalized=TRUE) %>% 
                     data.frame() %>% 
                     rownames_to_column(var = "gene")

tx2gene <- read.table(file = "data/tx2gene_grch38_ens94.txt", sep="\t", header = TRUE)

grch38annot <- tx2gene %>% 
  dplyr::select(ensgene, symbol) %>% 
  dplyr::distinct()

normalized_counts <- merge(normalized_counts, grch38annot, by.x="gene", by.y="ensgene")

normalized_counts <- normalized_counts %>%
  as_tibble()

#alternative----------------------------------------

#normalized_counts <- counts(dds, normalized=T) %>% 
#  data.frame() %>%
#  rownames_to_column(var="gene") %>%
#  as_tibble() %>%
#  left_join(grch38annot, by=c("gene" = "ensgene"))

grch38annot[grch38annot$symbol == "VPS35", "ensgene"]

# plotCounts(dds, gene="ENSG00000069329", intgroup="Genotype")
# #second-round---------
# d <- plotCounts(dds, gene="ENSG00000069329", intgroup="Genotype", returnData = TRUE)
# d %>% View()

d2 <- plotCounts(dds1, gene="ENSG00000069329", intgroup="Genotype", returnData = TRUE)
d2 %>% View()

d3 <- plotCounts(dds2, gene="ENSG00000069329", intgroup="Genotype", returnData = TRUE)
d3 %>% View()

# norm_sig1 <- normalized_counts[,c(1:13)] %>% 
#   filter(gene %in% sig1$gene)
# 
# norm_sig2 <- normalized_counts[,c(1,8:19)] %>% 
#   filter(gene %in% sig2$gene)

norm_sig3 <- normalized_counts[,c(1:13)] %>% 
  filter(gene %in% sig3$gene)

norm_sig4 <- normalized_counts[,c(1,8:19)] %>% 
  filter(gene %in% sig4$gene)

#setting the rownames to gene symbols----------------------------------
# norm_sig1 <- merge(norm_sig1, grch38annot, by.x="gene", by.y="ensgene")
# norm_sig1 <- norm_sig1 %>% column_to_rownames(var = "symbol")
# 
# norm_sig2 <- merge(norm_sig2, grch38annot, by.x="gene", by.y="ensgene")
# norm_sig2 <- norm_sig2 %>% column_to_rownames(var = "symbol")

norm_sig3 <- merge(norm_sig3, grch38annot, by.x="gene", by.y="ensgene")
norm_sig3 <- norm_sig3 %>% column_to_rownames(var = "symbol")

norm_sig4 <- merge(norm_sig4, grch38annot, by.x="gene", by.y="ensgene")
norm_sig4 <- norm_sig4 %>% column_to_rownames(var = "symbol")
#----------------------------------------------------------------------

#-------------------------HEATMAPS----------------------------------
heat_colors <- colorRampPalette(brewer.pal(6, "RdPu"))(100)

#----------------------------RES1 AND RES2 HEATMAPS - DON'T RUN-------------------------#
# pheatmap(norm_sig1[2:13], 
#          color = heat_colors, 
#          cluster_rows = T, 
#          show_rownames = F,
#          annotation = metadata2, 
#          border_color = NA, 
#          fontsize = 10, 
#          scale = "row", 
#          fontsize_row = 10, 
#          height = 20)
# 
# pheatmap(norm_sig2[2:13], 
#          color = heat_colors, 
#          cluster_rows = T, 
#          show_rownames = F,
#          annotation = metadata2, 
#          border_color = NA, 
#          fontsize = 10, 
#          scale = "row", 
#          fontsize_row = 10, 
#          height = 20)
#--------------------------------------------------------------------------------------#

#-----------------------VPS KO / WT heatmap-----------------------------------#
pheatmap(norm_sig3[2:13], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         annotation = metadata2, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 6, 
         height = 20,
         cutree_rows = 2)
#-----------------------VPS KO / VPS Rescue heatmap-----------------------------------#
pheatmap(norm_sig4[2:13], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         annotation = metadata2, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 6, 
         height = 20,
         cutree_rows = 2)
#-------------------------------------------------------------------------------------#
#------------------------------------RES1-(don't do this)------------------------------
# res1_tb <- res1_tb %>% 
#   mutate(threshold = padj < 0.05 & abs(log2FoldChange >=0.58))
# 
# ggplot(res1_tb) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
#   ggtitle("VPS depletion") +
#   xlab("log2 fold change") + 
#   ylab("-log10 adjusted p-value") +
#   #scale_y_continuous(limits = c(0,50)) +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25)))  
# 
# res1_tb <- bind_cols(res1_tb, symbol=grch38annot$symbol[match(res1_tb$gene, grch38annot$ensgene)])
# 
# res1_tb <- res1_tb %>% mutate(genelabels = "")
# 
# res1_tb <- res1_tb %>% arrange(padj)
# 
# res1_tb$genelabels[1:10] <- as.character(res1_tb$symbol[1:10])
# 
# ggplot(res1_tb, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(colour = threshold)) +
#   geom_text_repel(aes(label = genelabels)) +
#   ggtitle("VPS") +
#   xlab("log2 fold change") + 
#   ylab("-log10 adjusted p-value") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25))) 
# 
# #---------------------------RES2 (don't do this)---------------------------------#
# res2_tb <- res2_tb %>% 
#   mutate(threshold = padj < 0.05 & abs(log2FoldChange >=0.58))
# 
# ggplot(res2_tb) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
#   ggtitle("VPS depletion") +
#   xlab("log2 fold change") + 
#   ylab("-log10 adjusted p-value") +
#   #scale_y_continuous(limits = c(0,50)) +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25)))  
# 
# res2_tb <- bind_cols(res2_tb, symbol=grch38annot$symbol[match(res2_tb$gene, grch38annot$ensgene)])
# 
# res2_tb <- res2_tb %>% mutate(genelabels = "")
# 
# res2_tb <- res2_tb %>% arrange(padj)
# 
# res2_tb$genelabels[1:10] <- as.character(res2_tb$symbol[1:10])
# 
# #to make the volcano plot more comprehensive we remove VPS35 gene cause it has extremely low p-value and that is why it's fucking up the plot
# res2_tb <- res2_tb %>% filter(!symbol %in% c('VPS35'))
# 
# ggplot(res2_tb, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(colour = threshold)) +
#   geom_text_repel(aes(label = genelabels)) +
#   ggtitle("VPS KO vs VPS-GFP Rescue") +
#   xlab("log2 fold change") + 
#   ylab("-log10 adjusted p-value") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25))) 
# 
# res2_tb <- res2_tb %>% filter(!symbol %in% c('VPS35'))
# 
# #---------------------------------------------------------------------------------#

#-------------------------------------VPS KO / WT - Volcano prep and plot------------------#
res3_tb <- res3_tb %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >=0.58)

ggplot(res3_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("VPS depletion") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

res3_tb <- bind_cols(res3_tb, symbol=grch38annot$symbol[match(res3_tb$gene, grch38annot$ensgene)])

res3_tb <- res3_tb %>% mutate(genelabels = "")

res3_tb <- res3_tb %>% arrange(padj)

res3_tb$genelabels[1:18] <- as.character(res3_tb$symbol[1:18])

g <- ggplot(res3_tb, aes(x = log2FoldChange, y = -log10(padj), name = symbol)) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels), alpha = 0.5) +
  geom_vline(xintercept = c(0.58, -0.58), linetype = "longdash", alpha = 0.3) +
  geom_hline(yintercept = 0.05, linetype = "longdash", alpha = 0.3) +
  scale_color_manual(values = c("aquamarine4", "coral")) +
  ggtitle("VPS KO / WT") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.2), hjust = 0.5),
        axis.title = element_text(size = rel(1.1))) 

ggplotly(g)

##-------------------------------------VPS KO / VPS Rescue - Volcano prep and plot------------------#
res4_tb <- res4_tb %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >=0.58)

ggplot(res4_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("VPS depletion") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

res4_tb <- bind_cols(res4_tb, symbol=grch38annot$symbol[match(res4_tb$gene, grch38annot$ensgene)])

res4_tb <- res4_tb %>% mutate(genelabels = "")

res4_tb <- res4_tb %>% arrange(padj)

res4_tb$genelabels[1:18] <- as.character(res4_tb$symbol[1:18])

res4_tb <- res4_tb %>% filter(!symbol %in% c('VPS35'))

ggplot(res4_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels), alpha = 0.5) +
  geom_vline(xintercept = c(0.58, -0.58), linetype = "longdash", alpha = 0.3) +
  geom_hline(yintercept = 0.05, linetype = "longdash", alpha = 0.3) +
  scale_color_manual(values = c("aquamarine4", "coral")) +
  ggtitle("VPS KO / VPS-GFP Rescue") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.2), hjust = 0.5),
        axis.title = element_text(size = rel(1.1)))

#------------------------------------------------------------------------------------------------#
#biomaRt-example----------------------
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
listAttributes(ensembl)
listFilters(ensembl)
# 
# annotation <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 
#                                    'end_position', 'strand', 'gene_biotype', 
#                                    'external_gene_name', 'description'),
#       filters = 'ensembl_gene_id',
#       values = res4_tb$gene,
#       mart = ensembl)
# 
# annotated_res4 <- left_join(res4_tb, annotation, by=c('gene'='ensembl_gene_id'))
#-------------------------------

#--------------------------------GO enrichment-(VPS KO / WT)----------------------------------#
ent_gene <- getBM(attributes = c('entrezgene_id'),
                    filters = 'ensembl_gene_id',
                    values = sig3_up$gene,
                    mart = ensembl)
ent_gene <- ent_gene$entrezgene_id
ent_gene <- as.character(ent_gene)

ent_uni <- getBM(attributes = c('entrezgene_id'),
                  filters = 'ensembl_gene_id',
                  values = res3_tb$gene,
                  mart = ensembl)
ent_uni <- ent_uni$entrezgene_id 
ent_uni <- as.character(ent_uni)

ego_up <- enrichGO(gene = ent_gene,
                OrgDb = org.Hs.eg.db,
                ont = 'BP',
                universe = ent_uni,
                readable = TRUE)
#------------------------------------------------------------------------------------#
#--------------------------------GO depletion-(VPS KO / WT)----------------------------------#
ent_gene2 <- getBM(attributes = c('entrezgene_id'),
                  filters = 'ensembl_gene_id',
                  values = sig3_down$gene,
                  mart = ensembl)
ent_gene2 <- ent_gene2$entrezgene_id
ent_gene2 <- as.character(ent_gene2)

ent_uni2 <- getBM(attributes = c('entrezgene_id'),
                 filters = 'ensembl_gene_id',
                 values = res3_tb$gene,
                 mart = ensembl)
ent_uni2 <- ent_uni2$entrezgene_id 
ent_uni2 <- as.character(ent_uni2)

ego_down <- enrichGO(gene = ent_gene2,
                OrgDb = org.Hs.eg.db,
                ont = 'BP',
                universe = ent_uni2,
                readable = TRUE)
#VISUALIZAING
clusterProfiler::cnetplot(ego_up, colorEdge = TRUE)
dotplot(ego_up)
goplot(ego_up)
clusterProfiler::cnetplot(ego_down)
dotplot(ego_down)
goplot(ego_down)
#----------------------------------------KEGG - VPS KO / WT-----------------------------------#
ekeg <- enrichKEGG(gene = ent_gene,
                 organism = "hsa",
                 universe = ent_uni)
clusterProfiler::cnetplot(ekeg)

#---------------------------------------------------------------------------------------------#
#--------------------------------GO enrichment-(VPS KO / VPS Rescue)----------------------------------#

ent_geneB <- getBM(attributes = c('entrezgene_id'),
                  filters = 'ensembl_gene_id',
                  values = sig4_up$gene,
                  mart = ensembl)
ent_geneB <- ent_geneB$entrezgene_id
ent_geneB <- as.character(ent_geneB)

ent_uniB <- getBM(attributes = c('entrezgene_id'),
                 filters = 'ensembl_gene_id',
                 values = res4_tb$gene,
                 mart = ensembl)
ent_uniB <- ent_uniB$entrezgene_id 
ent_uniB <- as.character(ent_uniB)

ego_upB <- enrichGO(gene = ent_geneB,
                   OrgDb = org.Hs.eg.db,
                   ont = 'BP',
                   universe = ent_uniB,
                   pvalueCutoff = 0.1,
                   readable = TRUE)
#--------------------------------------------------------------------------------------------#
#--------------------------------GO depletion-(VPS KO / VPS Rescue)----------------------------------#

ent_geneB2 <- getBM(attributes = c('entrezgene_id'),
                   filters = 'ensembl_gene_id',
                   values = sig4_down$gene,
                   mart = ensembl)
ent_geneB2 <- ent_geneB2$entrezgene_id
ent_geneB2 <- as.character(ent_geneB2)

ent_uniB2 <- getBM(attributes = c('entrezgene_id'),
                  filters = 'ensembl_gene_id',
                  values = res4_tb$gene,
                  mart = ensembl)
ent_uniB2 <- ent_uniB2$entrezgene_id 
ent_uniB2 <- as.character(ent_uniB2)

ego_downB <- enrichGO(gene = ent_geneB2,
                     OrgDb = org.Hs.eg.db,
                     ont = 'BP',
                     universe = ent_uniB2,
                     pvalueCutoff = 0.2,
                     readable = TRUE)

#VISUALIZING
clusterProfiler::cnetplot(ego_upB, colorEdge = TRUE)
dotplot(ego_upB)
goplot(ego_upB)
clusterProfiler::cnetplot(ego_downB)
dotplot(ego_downB)
goplot(ego_downB)
#----------------------------------------KEGG - VPS KO / VPS Rescue-----------------------------------#

ekegB <- enrichKEGG(gene = ent_geneB,
                   organism = "hsa",
                   universe = ent_uniB)
clusterProfiler::cnetplot(ekegB)

#----------------------------------CORRELATION TESTS--(KO/WT vs KO/Rescue)-------------------------------------------------#

sig3_gene <- sig3$gene
sig4_gene <- sig4$gene

common_gene <- intersect(sig3_gene, sig4_gene)

sig3_common <- sig3[sig3$gene %in% common_gene, c('gene', 'log2FoldChange')]
sig4_common <- sig4[sig4$gene %in% common_gene, c('gene', 'log2FoldChange')]

common_fold_changes <- left_join(sig3_common, sig4_common, by = c('gene' = 'gene'))

common_fold_changes <- common_fold_changes %>% 
  mutate(upregulated = log2FoldChange.x > 0) 

ggplot(data = common_fold_changes, aes(x=log2FoldChange.x, y=log2FoldChange.y, color = upregulated)) +
  geom_point() +
  scale_color_manual(values=c("coral","aquamarine4")) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme_bw()

cor.test(x = common_fold_changes$log2FoldChange.x, 
         y = common_fold_changes$log2FoldChange.y,
         method = 'spearman')

cor.test(x = common_fold_changes$log2FoldChange.x, 
         y = common_fold_changes$log2FoldChange.y,
         method = 'pearson')
#------------------------------------VENN DIAGRAM---(KO/WT vs KO/Rescue)---------------------------------------------#

sig3_only <- setdiff(sig3_gene, sig4_gene)
sig4_only <- setdiff(sig4_gene, sig3_gene)

plot.new()
draw.pairwise.venn(area1 = length(sig3_gene),
                   area2 = length(sig4_gene),
                   cross.area = length(common_gene),
                   scaled = TRUE, fill= c('red', 'purple'))
