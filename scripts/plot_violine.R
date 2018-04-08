.libPaths(rev(.libPaths()))
library(plyr)
library(dplyr) # Dataframe manipulation
library(Matrix) # Sparse matrices
library(stringr)
library(RColorBrewer)
#library(matrixStats)
# library(Hmisc) # for cut2 function
library(devtools)
library(Seurat)
library(plotly)

# importing UMI
# umi_matrix             <- read.csv(file.path(path,'summary/umi_expression_matrix.tsv'), sep='\t', header = TRUE, row.names = 1, check.names = FALSE)
headm <- function(x, n=6) x[1:6,1:6]
count_matrix             <- read.csv(snakemake@input$counts, sep='\t', header = TRUE, row.names = 1, check.names = FALSE)
# importing counts
# summary/counts_expression_matrix.tsv
umi_matrix           <- read.csv(snakemake@input$UMIs, sep='\t', header = TRUE, row.names = 1, check.names = FALSE)
# colnames(count_matrix) <- colnames(umi_matrix)
design                 <- read.csv(snakemake@input$design, stringsAsFactors = TRUE, header = TRUE, row.names = NULL)

metaData <- data.frame(cellNames = colnames(umi_matrix)) %>%
  mutate(samples = factor(str_replace(cellNames,"_[^_]*$",""))) %>%
  mutate(barcode = factor(str_replace(cellNames,".+_",""))) %>%
  left_join(design, by = "samples")
rownames(metaData) <- metaData$cellNames

# setting is.expr = -1 to avoid filtering whilst creating
# myumi <- CreateSeuratObject(raw.data = umi_matrix, meta.data = metaData, is.expr = -1)
myumi <- CreateSeuratObject(raw.data = umi_matrix, meta.data = metaData)
myumi <- SetAllIdent(object = myumi, id = "samples")
# relabel cell idenity (https://github.com/satijalab/seurat/issues/380)
# myumi@ident  <- factor(stringr::str_replace(colnames(myumi@data),"_[^_]+$",""))
# myumi@meta.data$myident  <- (stringr::str_replace(colnames(myumi@data),"_[^_]+$",""))

mycount <- CreateSeuratObject(raw.data = count_matrix, meta.data = metaData)
mycount <- SetAllIdent(object = mycount, id = "samples")
# turn off filtering

# note, the @meta.data slot contains usefull summary stuff 
# head(mycount@meta.data,2)
#                              nGene nUMI expected_cells read_length      barcode
# dropseqLib1_ACTAACATTATT    15   33            400         100 ACTAACATTATT
# dropseqLib1_GAGTCTGAGGCG     5    9            400         100 GAGTCTGAGGCG
#                                       origin      origin
# dropseqLib1_ACTAACATTATT dropseqLib1 dropseqLib1
# dropseqLib1_GAGTCTGAGGCG dropseqLib1 dropseqLib1
meta.data         <- myumi@meta.data
meta.data$nCounts <- mycount@meta.data$nUMI
# save.image(file="Rworkspace_violine.rdata")
setwd("~/analysis/dropseq/data/10x")
# load("Rworkspace_violine.rdata")
load("Seurat_objects.rdata")

gg <- ggplot(meta.data, aes(x=nUMI,y=nCounts,color=orig.ident)) +
  #   coord_trans(y="log10",x = "log10") +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
plotname <- "umi_vs_counts"
p <- ggplotly(gg)
htmlwidgets::saveWidget(p,file.path(paste0(plotname,".html")))

# how about unaligned reads/UMI?
# Note(Seb): raw.data is actually filtered data i.e. nr of genes likely to be smaller than input data!
mito.gene.names  <- grep("^mt-", rownames(myumi@raw.data), value=TRUE)
# mito.gene.names2  <- subset(mito.gene.names, mito.gene.names %in% rownames(myumi@raw.data))
sribo.gene.names <- grep("^Rps", rownames(myumi@raw.data), value=TRUE)
lribo.gene.names <- grep("^Rpl", rownames(myumi@raw.data), value=TRUE)

col.total.umi            <- Matrix::colSums(myumi@raw.data)
col.total.count          <- Matrix::colSums(mycount@raw.data)
meta.data$col.total.umi   <- col.total.umi
meta.data$col.total.count <- col.total.count
# mito.percent.counts <- mito.percent.counts
# sribo.pct <- sribo.pct
# lribo.pct <- lribo.pct
# ribo_tot <- ribo_tot

myumi.top_50   <- apply(myumi@raw.data, 2, function(x) sum(x[order(x, decreasing = TRUE)][1:50])/sum(x))
mycount.top_50 <- apply(mycount@raw.data, 2, function(x) sum(x[order(x, decreasing = TRUE)][1:50])/sum(x))

myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[sribo.gene.names, ])/col.total.umi, "pct.sribo")
myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[lribo.gene.names, ])/col.total.umi, "pct.lribo")
myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[unique(c(sribo.gene.names, lribo.gene.names)), ])/col.total.umi, "pct.Ribo")
myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[mito.gene.names, ])/col.total.umi, "pct.mito")
myumi <- AddMetaData(myumi, myumi.top_50, 'top50')
tmp <- myumi@meta.data$nUMI/myumi@meta.data$nGene
names(tmp) <- rownames(myumi@meta.data)
myumi <- AddMetaData(myumi, tmp, 'umi.per.gene')

tmp=(myumi@meta.data[,"nUMI",drop=F]/myumi@meta.data$nGene)
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[sribo.gene.names, ])/col.total.count, "pct.sribo")
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[lribo.gene.names, ])/col.total.count, "pct.lribo")
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[unique(c(sribo.gene.names, lribo.gene.names)), ])/col.total.count, "pct.Ribo")
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[mito.gene.names, ])/col.total.count, "pct.mito")
mycount <- AddMetaData(mycount, mycount.top_50, 'top50')
tmp <- mycount@meta.data$nUMI/mycount@meta.data$nGene
names(tmp) <- rownames(mycount@meta.data)
mycount <- AddMetaData(mycount, tmp, 'umi.per.gene')
# mycount@meta.data$count.per.gene <- mycount@meta.data$nUMI/mycount@meta.data$nGene

# mytheme <- theme_bw(base_size = 9) +
mytheme <- theme_bw() +
  theme(legend.position = "right",
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 300, hjust = 0))
theme_set(mytheme)

gg <- VlnPlot(myumi,c("nUMI", "nGene", "top50", 'umi.per.gene','pct.Ribo', "pct.mito"), x.lab.rot = TRUE, do.return = TRUE)
# ggsave(gg,file=file.path("violinplots_comparison_UMI.pdf"),width=18,height=18)
ggsave(gg, file=snakemake@output$pdf, width = 18, height = 18)
# gg <- VlnPlot(mycount,c("nUMI", "nGene", "top50", 'count.per.gene','pct.Ribo', "pct.mito"), x.lab.rot = TRUE, do.return = TRUE)
# ggsave(gg,file=file.path("violinplots_comparison_count.pdf"),width=18,height=18)

# gg <- GenePlot(object = myumi, gene1 = "nUMI", gene2 = "nGene")
# ggsave(gg,file=file.path("violinplots_comparison.pdf"),width=18,height=18)

gg <- ggplot(myumi@meta.data, aes(x = nGene, y = nUMI, color=orig.ident)) +
  geom_point(size=.5) +
  geom_smooth() +
  labs(title = "Genes (pooled mouse and human set) vs UMIs for each bead",
       x = "Number of UMIs per Bead",
       y = "Number of Genes per Bead")
plotname <- "umi_vs_gene"
p <- ggplotly(gg)
htmlwidgets::saveWidget(p,file.path(paste0(plotname,".html")))
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)

gg <- ggplot(myumi@meta.data, aes(x = nGene, y = nUMI, color=orig.ident)) +
  geom_point(size=.5) +
  geom_smooth() +
  xlim(0,4e4) +
  ylim(0,10e3) +
  labs(title = "Genes (pooled mouse and human set) vs UMIs for each bead",
       x = "Number of UMIs per Bead",
       y = "Number of Genes per Bead")
plotname <- "umi_vs_gene_zoom"
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)
gg <- ggplot(myumi@meta.data, aes(x = nGene, y = nUMI, color=orig.ident)) +
  geom_point(size=.5) +
  geom_smooth() +
  coord_trans(x = "log10") +
#   scale_y_continuous(trans='log10') +
  labs(title = "Genes (pooled mouse and human set) vs UMIs for each bead",
       x = "Number of UMIs per Bead",
       y = "Number of Genes per Bead")
plotname <- "umi_vs_gene_log"
p <- ggplotly(gg)
htmlwidgets::saveWidget(p,file.path(paste0(plotname,".html")))
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)
gg <- ggplot(mycount@meta.data, aes(x = nGene, y = nUMI, color=orig.ident)) +
  geom_smooth() +
  geom_point(size=.5) +
  labs(title = "Genes (pooled mouse and human set) vs read counts for each bead",
       x = "Number of Reads per Bead",
       y = "Number of Genes per Bead")
plotname <- "count_vs_gene"
p <- ggplotly(gg)
htmlwidgets::saveWidget(p,file.path(paste0(plotname,".html")))
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)
gg <- ggplot(mycount@meta.data, aes(x = nGene, y = nUMI, color=orig.ident)) +
  geom_smooth() +
  geom_point(size=.5) +
  xlim(0,1e5) +
  ylim(0,10e3) +
  labs(title = "Genes (pooled mouse and human set) vs read counts for each bead",
       x = "Number of Reads per Bead",
       y = "Number of Genes per Bead")
plotname <- "count_vs_gene_zoom"
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)
gg <- ggplot(mycount@meta.data, aes(x = nGene, y = nUMI, color=orig.ident)) +
  geom_smooth() +
  geom_point(size=.5) +
  coord_trans(x = "log10") +
#  scale_y_continuous(trans='log10') +
  labs(title = "Genes (pooled mouse and human set) vs read counts for each bead",
       x = "Number of Reads per Bead",
       y = "Number of Genes per Bead")
plotname <- "count_vs_gene_log"
p <- ggplotly(gg)
htmlwidgets::saveWidget(p,file.path(paste0(plotname,".html")))
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)

save.image(file="Seurat_objects.rdata")
