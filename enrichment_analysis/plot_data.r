setwd("~/Documents/PhD/rotations/KoremRotation/enrichment_analysis/")
library(ggplot2)
library(grid)
library(gridExtra)


data_CP001509    <- read.table("CP001509.offset.gff3")
data_CP001637    <- read.table("CP001637.offset.gff3")
data_NC_009800   <- read.table("NC_009800.offset.gff3")
data_NC_009801   <- read.table("NC_009801.offset.gff3")
data_NC_010468   <- read.table("NC_010468.offset.gff3")
data_NC_010498   <- read.table("NC_010498.offset.gff3")
data_NC_011353   <- read.table("NC_011353.offset.gff3")
data_NC_012892   <- read.table("NC_012892.offset.gff3")


CP001509__plot = ggplot(data_CP001509 ,aes(x=V10, y=V11)) + geom_boxplot(aes(color=V10)) + ggtitle("CP001509") 
CP001637__plot = ggplot(data_CP001637 ,aes(x=V10, y=V11)) + geom_boxplot(aes(color=V10)) + ggtitle("CP001637") 
NC_009800_plot = ggplot(data_NC_009800,aes(x=V10, y=V11)) + geom_boxplot(aes(color=V10)) + ggtitle("NC_009800")
NC_009801_plot = ggplot(data_NC_009801,aes(x=V10, y=V11)) + geom_boxplot(aes(color=V10)) + ggtitle("NC_009801")
NC_010468_plot = ggplot(data_NC_010468,aes(x=V10, y=V11)) + geom_boxplot(aes(color=V10)) + ggtitle("NC_010468")
NC_010498_plot = ggplot(data_NC_010498,aes(x=V10, y=V11)) + geom_boxplot(aes(color=V10)) + ggtitle("NC_010498")
NC_011353_plot = ggplot(data_NC_011353,aes(x=V10, y=V11)) + geom_boxplot(aes(color=V10)) + ggtitle("NC_011353")
NC_012892_plot = ggplot(data_NC_012892,aes(x=V10, y=V11)) + geom_boxplot(aes(color=V10)) + ggtitle("NC_012892")
grid.arrange(
  CP001509__plot,
  CP001637__plot,
  NC_009800_plot,
  NC_009801_plot,
  NC_010468_plot,
  NC_010498_plot,
  NC_011353_plot,
  NC_012892_plot,
  ncol=3, nrow=3,
  top = textGrob("Distribution of genes lethal at high copy number."),
  left = textGrob("Position with respect to Ori (bp)", rot=90),
  bottom= textGrob("PanDaTox (TRUE) non-PanDaTox(FALSE)"))


# perform wilcoxon rank sum test to see if PanDaTox genes are significantly
# greater than non-pandatox genes
sink('wilcoxon_results.txt')
wilcox.test(data_CP001509$V11  ~ data_CP001509$V10, alternative="less")
wilcox.test(data_CP001637$V11  ~ data_CP001637$V10, alternative="less")
wilcox.test(data_NC_009800$V11 ~ data_NC_009800$V10, alternative="less")
wilcox.test(data_NC_009801$V11 ~ data_NC_009801$V10, alternative="less")
wilcox.test(data_NC_010468$V11 ~ data_NC_010468$V10, alternative="less")
wilcox.test(data_NC_010498$V11 ~ data_NC_010498$V10, alternative="less")
wilcox.test(data_NC_011353$V11 ~ data_NC_011353$V10, alternative="less")
wilcox.test(data_NC_012892$V11 ~ data_NC_012892$V10, alternative="less")
sink()
