# modified by Sebastian Mueller (sebm_at_posteo.de)
# from rules/map.smk:
#rule plot_knee_plot:
#	input: 'logs/{sample}_hist_out_cell.txt'
#	output:
#		pdf = 'plots/{sample}_knee_plot.pdf',
#		png = 'plots/png/{sample}_knee_plot.png'
#	params: 
#		cells = lambda wildcards: samples.loc[wildcards.sample,'expected_cells']
#	script:
#		'../scripts/plot_knee_plot.R'

# input format of hist_out_cell.txt (each line):
# count \tab cell-barcode(12bp)
# counts are not filtered for PCR duplicates
library(ggplot2)
library(plyr)
library(gridExtra)
library(grid)
# potential usefull to change lag of diff calculation
mylag <- 1
data <- read.table(file = snakemake@input[[1]][1], header=FALSE, stringsAsFactors=FALSE)
# data <- read.table(file = file, header=FALSE, stringsAsFactors=FALSE)
# head(data,3)
#        V1           V2
# 1 1145137 CCCTTCGTCTGC
# 2 1039974 ATAGTTTTTTAA
# 3  912199 GCATGAAACTTC

# borrowed from https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

reads        <- data$V1
barcodes     <- data$V2
total_reads  <- sum(reads)
reads_cumsum <- cumsum(reads)
# 1st dervivative (diff) needs also to be padded to keep same vector length 
reads_diff   <- c(diff(reads,lag=mylag,differences = 1),rep(0,mylag))

# 2nd derivative: twice as much padding:
reads_diff_diff   <- c(diff(reads,lag=mylag,differences = 2),rep(0,mylag*2))
reads_cumsum_perc <- (reads_cumsum/total_reads)
x                 <- 1:length(reads_cumsum_perc)
plot_data         <- data.frame(rank = x,cum_sum=reads_cumsum_perc,read_count=reads, Barcode=data$V2, diff=reads_diff,diffdiff=reads_diff_diff)

# head(plot_data,6)
#   rank     cum_sum read_count      Barcode    diff diffdiff
# 1    1 0.007192916    1145137 CCCTTCGTCTGC -105163   -22612
# 2    2 0.013725274    1039974 ATAGTTTTTTAA -127775    14025
# 3    3 0.019455043     912199 GCATGAAACTTC -113750    61486
# 4    4 0.024470318     798449 GTGTGGGTCTCT  -52264    34985
# 5    5 0.029157308     746185 CGTACTGACTAC  -17279   -60441
# 6    6 0.033735764     728906 GTTCGTCCCGCC  -77720    69104
mystats <- paste0("| Nr barcodes total: ", length(barcodes), ' \n ',
                 "Nr barcodes for 50% reads: ", which.min(reads_cumsum_perc<0.50)," | ",
                 "Nr barcodes for 95% reads: ", which.min(reads_cumsum_perc<0.95)," | ",
                 "Nr barcodes for 99% reads: ", which.min(reads_cumsum_perc<0.99)
                 )

# x_scale        <- which(reads_cumsum_perc>0.99)[1]
x_scale        <- snakemake@params$cells * 4
plot_data_head <- head(plot_data, x_scale)
# plot_data_head$reads_diff_smooth <- predict(loess(diff~rank,data=plot_data_head))
plot_data_head$reads_diffdiff_smooth <- predict(loess(diffdiff~rank,span=0.2,data=plot_data_head))

## Finding knee in knee-plot:
# Best approach so far is to calclate the 2nd derivative of the read counts per STAMP (reads_diff_diff), smooth it (loess) and find it's maxima. Since ther can be several maxima, they are just ploted and it's up to the user to visually asses and decide.
# finding local maxima in 2nd derivative:
#https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
loc.max <- localMaxima(plot_data_head$reads_diffdiff_smooth)
# local maxima are then colored green as lines in plots 
plot_data_head_sub <- plot_data_head[loc.max,]
# which.peaks(plot_data_head$reads_diff_smooth2)

knee_plot <- ggplot(plot_data_head, aes(x=rank, y=cum_sum)) +
  xlim(0,x_scale) +
  ylim(0,1) +
  geom_text(data=plot_data_head_sub,aes(label = round(cum_sum,2)),nudge_y=-0.05, vjust = "inward", hjust = "inward") +
	geom_vline(xintercept=snakemake@params$cells, linetype="dashed", color = "red") +
	#   geom_vline(xintercept=100, linetype="dashed", color = "red") +
  geom_vline(xintercept=loc.max, col="lightgreen") +
  geom_point(size = 0.1)  +
	ggtitle(paste0(snakemake@wildcards$sample, '\nTotal reads: ', prettyNum(total_reads), mystats)) +
  theme(plot.title = element_text(size=10)) +
  labs(x='STAMPS', y='Cumulative fraction of reads')
read_count_plot <- ggplot(plot_data_head, aes(x=rank, y=read_count)) +
  geom_text(data=plot_data_head_sub,aes(label = read_count),nudge_y=-0.05, vjust = "inward", hjust = "inward", check_overlap = TRUE) +
  #   geom_smooth() +
  xlim(0,x_scale) +
  #   ylim(0,1000) +
	geom_vline(xintercept=snakemake@params$cells, linetype="dashed", color = "red") +
	#   geom_vline(xintercept=100, linetype="dashed", color = "red") +
  geom_vline(xintercept=loc.max, col="lightgreen") +
  geom_point(size = 0.1)  +
  theme(plot.title = element_text(size=10)) +
  labs(x='STAMPS', y='Read counts per STAMP')
diff_plot <- ggplot(plot_data_head, aes(x=rank, y=diff)) +
  geom_text(data=plot_data_head_sub,aes(label = diff),nudge_y=-0.05, vjust = "inward", hjust = "inward", check_overlap = TRUE) +
  #   geom_smooth() +
  xlim(0,x_scale) +
  #   ylim(0,1000) +
	geom_vline(xintercept=snakemake@params$cells, linetype="dashed", color = "red") +
	#   geom_vline(xintercept=100, linetype="dashed", color = "red") +
  geom_vline(xintercept=loc.max, col="lightgreen") +
  geom_point(size = 0.1)  +
  theme(plot.title = element_text(size=10)) +
  #1st derivative read count diff to next STAMP")
  labs(x='STAMPS', y="1st derivative of read counts")
diff_diff_plot <- ggplot(plot_data_head, aes(x=rank, y=reads_diffdiff_smooth)) +
  geom_text(data=plot_data_head_sub,aes(label = rank),nudge_y=-1, vjust = "inward", hjust = "inward", check_overlap = TRUE) +
  xlim(0,x_scale) +
  ylim(-30,30) +
	geom_vline(xintercept=snakemake@params$cells, linetype="dashed", color = "red") +
	#   geom_vline(xintercept=100, linetype="dashed", color = "red") +
  geom_vline(xintercept=loc.max, col="lightgreen") +
  geom_point(size = 0.1)  +
  theme(plot.title = element_text(size=10)) +
  labs(x='STAMPS', y='2nd derivative of read counts')
	#   scale_y_continuous(position = "right")
gp1 <- ggplotGrob(knee_plot)
gp2 <- ggplotGrob(read_count_plot)
gp3 <- ggplotGrob(diff_plot)
gp4 <- ggplotGrob(diff_diff_plot)
grid::grid.newpage()
gg=grid::grid.draw(rbind(gp1, gp2, gp3, gp4, size = "last"))
gg=gridExtra::arrangeGrob(rbind(gp1, gp2, gp3, gp4, size = "last"))

# if barcode.csv is present in base directory, only use barcodes in there (rule plot_knee_plot_whitelist in map.smk)
if(!is.null(snakemake@input$barcodes))
{
  barcodes_whitelist  <-  read.csv(snakemake@input$barcodes, header=FALSE, stringsAsFactors=FALSE)
  knee_plot  <-  knee_plot +
    geom_point(data = plot_data[plot_data$Barcode %in% barcodes_whitelist$V1,],
               aes(x=rank, y=cum_sum, color='Whitelisted'), size=0.1) +
    scale_color_manual(values=c('Whitelisted'='green'))
}
ggsave(gg, file=snakemake@output$pdf, width = 9, height = 11)
# ggsave(knee_plot, file=snakemake@output$pdf, width = 4, height = 3)
ggsave(gg, file=snakemake@output$png, width = 9, height = 11)

# suggestions:
# plot deduplicated reads instead
# note: UMI count (such as in cell ranger) in might not be sufficient to mark PCR duplicates, the associate RNA sequence has to be duplicated as well
# color points according to some quality value, i.e. mapping rate etc.
