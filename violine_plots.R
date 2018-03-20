library(plyr)
library(dplyr) # Dataframe manipulation
library(Matrix) # Sparse matrices
library(stringr)
library(RColorBrewer)
#library(matrixStats)
library(Hmisc) # for cut2 function
library(devtools)
library(Seurat)
library(plotly)
# known bugs:
# problems if sample names have underscores in its name

umi_matrix             <- read.csv(file.path(path,'summary/umi_expression_matrix2.tsv'), sep='\t', header = TRUE, row.names = 1, check.names = FALSE)
count_matrix           <- read.csv(file.path(path,'summary/counts_expression_matrix.tsv'), sep='\t', header = TRUE, row.names = 1, check.names = FALSE)
colnames(count_matrix) <- colnames(umi_matrix)
design                 <- read.csv(file.path(path,'samples2.csv'), stringsAsFactors = TRUE, header = TRUE, row.names = 1)

CreatePdata <- function(pdata, cell.names) {
  for (i in 1:ncol(pdata)){
    pdata[,i] <- factor(pdata[,i])
  }
  fullMeta <- data.frame(matrix(ncol = ncol(pdata), nrow = length(cell.names)))
  colnames(fullMeta) <- colnames(pdata)
  #   samples <- factor(str_split(cell.names, '_', simplify = TRUE)[,1])
  samples <- factor(str_replace(cell.names,"_[^_]*$",""))
  for (i in 1:ncol(pdata)){
    column_values <- mapvalues(samples,
                              from = levels(samples), 
                              to = pdata[levels(samples),i])
    column_values <- mapvalues(column_values,
                              from=levels(column_values),
                              to=levels(pdata[,i]))
    fullMeta[,i] <- column_values
  }
  rownames(fullMeta) <- cell.names
  return(fullMeta)
}
fullMeta         <- CreatePdata(design, cell.names = colnames(umi_matrix))
fullMeta$barcode <- str_replace(rownames(fullMeta),".+_","")
fullMeta$origin  <- str_replace(rownames(fullMeta),"_.+","")

# myumi <- CreateSeuratObject(raw.data = umi_matrix[,-test], meta.data = fullMeta)
myumi <- CreateSeuratObject(raw.data = umi_matrix, meta.data = fullMeta)
mycount <- CreateSeuratObject(raw.data = count_matrix, meta.data = fullMeta)
# note, the @meta.data slot contains usefull summary stuff 
# head(mycount@meta.data,2)
#                              nGene nUMI expected_cells read_length      barcode
# dropseqLib1_ACTAACATTATT    15   33            400         100 ACTAACATTATT
# dropseqLib1_GAGTCTGAGGCG     5    9            400         100 GAGTCTGAGGCG
#                                       origin      orig.ident
# dropseqLib1_ACTAACATTATT dropseqLib1 dropseqLib1
# dropseqLib1_GAGTCTGAGGCG dropseqLib1 dropseqLib1
head(mycount@meta.data,2)
meta.data         <- myumi@meta.data
meta.data$nCounts <- mycount@meta.data$nUMI

gg <- ggplot(meta.data, aes(x=nUMI,y=nCounts,color=orig.ident)) +
  #   coord_trans(y="log10",x = "log10") +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
plotname <- "umi_vs_counts"
p <- ggplotly(gg)
htmlwidgets::saveWidget(p,file.path(paste0(plotname,".html")))

# ribo_genes <- read.table('2017/BigExp/mart_export.txt', sep = '\t', header = TRUE)
# all_ribo.genes <- unique(ribo_genes$Gene.name)
# final_list=NULL
# for (i in unlist(all_ribo.genes)){
# print(i)
#   final_list <- c(final_list,(grep(pattern = paste0(i,".*"), x = rownames(myumi@raw.data), value = TRUE)))
# }
# print(length(final_list))

# how about unaligned reads/UMI?
mito.gene.names  <- grep("^mt-", rownames(myumi@raw.data), value=TRUE)
mito.gene.names  <- subset(mito.gene.names, mito.gene.names %in% rownames(myumi@raw.data))
sribo.gene.names <- c(grep("^Rps", rownames(myumi@raw.data), value=TRUE))
lribo.gene.names <- c(grep("^Rpl", rownames(myumi@raw.data), value=TRUE))

col.total.umi            <- Matrix::colSums(myumi@raw.data)
col.total.count          <- Matrix::colSums(mycount@raw.data)
fullMeta$col.total.umi   <- col.total.umi
fullMeta$col.total.count <- col.total.count
# mito.percent.counts <- mito.percent.counts
# sribo.pct <- sribo.pct
# lribo.pct <- lribo.pct
# ribo_tot <- ribo_tot

myumi.top_50   <- apply(myumi@raw.data, 2, function(x) sum(x[order(x, decreasing = TRUE)][1:50])/sum(x))
mycount.top_50 <- apply(mycount@raw.data, 2, function(x) sum(x[order(x, decreasing = TRUE)][1:50])/sum(x))

myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[sribo.gene.names, ])/col.total.umi, "pct.sribo")
myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[lribo.gene.names, ])/col.total.umi, "pct.lribo")
# myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[unique(c(final_list,sribo.gene.names, lribo.gene.names)), ])/col.total.umi, "pct.Ribo")
myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[mito.gene.names, ])/col.total.umi, "pct.mito")
myumi <- AddMetaData(myumi, myumi.top_50, 'top50')
myumi@meta.data$umi.per.gene <- myumi@meta.data$nUMI/myumi@meta.data$nGene

mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[sribo.gene.names, ])/col.total.count, "pct.sribo")
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[lribo.gene.names, ])/col.total.count, "pct.lribo")
# mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[unique(c(final_list,sribo.gene.names, lribo.gene.names)), ])/col.total.count, "pct.Ribo")
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[mito.gene.names, ])/col.total.count, "pct.mito")
mycount <- AddMetaData(mycount, mycount.top_50, 'top50')
mycount@meta.data$count.per.gene=mycount@meta.data$nUMI/mycount@meta.data$nGene


mytheme <- theme_bw(base_size = 9) + 
  theme(legend.position = "right",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 300, hjust = 0))
theme_set(mytheme)

gg <- VlnPlot(myumi,c("nUMI", "nGene", "top50", 'umi.per.gene','pct.Ribo', "pct.mito"), x.lab.rot = TRUE, do.return = TRUE)
ggsave(gg,file=file.path("violinplots_comparison_UMI.pdf"),width=18,height=18)
# gg <- VlnPlot(mycount,c("nUMI", "nGene", "top50", 'count.per.gene','pct.Ribo', "pct.mito"), x.lab.rot = TRUE, do.return = TRUE)
# ggsave(gg,file=file.path("violinplots_comparison_count.pdf"),width=18,height=18)

GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = myumi, gene1 = "nUMI", gene2 = "nGene")
ggsave(gg,file=file.path("violinplots_comparison.pdf"),width=18,height=18)

gg <- ggplot(myumi@meta.data, aes(x = nUMI, y = nGene, color=orig.ident)) +
	geom_point(size=.5) +
	geom_smooth() +
	labs(title = "Genes (pooled mouse and human set) vs UMIs for each bead",
			 x = "Number of UMIs per Bead",
			 y = "Number of Genes per Bead")
plotname <- "umi_vs_gene"
p <- ggplotly(gg)
htmlwidgets::saveWidget(p,file.path(paste0(plotname,".html")))
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)
gg <- ggplot(myumi@meta.data, aes(x = nUMI, y = nGene, color=orig.ident)) +
	geom_point(size=.5) +
	geom_smooth() +
	xlim(0,4e4) +
	ylim(0,10e3) +
	labs(title = "Genes (pooled mouse and human set) vs UMIs for each bead",
			 x = "Number of UMIs per Bead",
			 y = "Number of Genes per Bead")
plotname <- "umi_vs_gene_zoom"
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)
gg <- ggplot(myumi@meta.data, aes(x = nUMI, y = nGene, color=orig.ident)) +
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
gg <- ggplot(mycount@meta.data, aes(x = nUMI, y = nGene, color=orig.ident)) +
	geom_smooth() +
	geom_point(size=.5) +
	labs(title = "Genes (pooled mouse and human set) vs read counts for each bead",
			 x = "Number of Reads per Bead",
			 y = "Number of Genes per Bead")
plotname <- "count_vs_gene"
p <- ggplotly(gg)
htmlwidgets::saveWidget(p,file.path(paste0(plotname,".html")))
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)
gg <- ggplot(mycount@meta.data, aes(x = nUMI, y = nGene, color=orig.ident)) +
	geom_smooth() +
	geom_point(size=.5) +
	xlim(0,1e5) +
	ylim(0,10e3) +
	labs(title = "Genes (pooled mouse and human set) vs read counts for each bead",
			 x = "Number of Reads per Bead",
			 y = "Number of Genes per Bead")
plotname <- "count_vs_gene_zoom"
ggsave(gg,file=file.path(paste0(plotname,".pdf")),width=12,height=7)
gg <- ggplot(mycount@meta.data, aes(x = nUMI, y = nGene, color=orig.ident)) +
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

save.image(file="Rworkspace.rdata")
