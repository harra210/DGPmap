#!/usr/bin/env Rscript
## This script is modified from the original creator of this script Dr. Reuben Buckley, Ph.D.. ##

## Packages required to run script ##
library(optparse)
library(matrixStats)
library(zoo)
library(UpSetR)
##
#
rm(list = ls()) # this clears out the environment and ensure a clean slate
options(stringsAsFactors = FALSE)
#
option_list = list(
  make_option(c("-i","--input"), type="character", default=NULL, help="Parent directory containing QC files", metavar="character"),
      make_option(c("-d","--dataset"), type="character", default=NULL, help="Name of the dataset being analyzed (eg. 3176g_AutoXPAR)", metavar="character"),
        make_option(c("-f","--filter"), type="character", default=NULL, help="File containing sample names to be filtered out of analysis, if there are none supply an empty text file.", metavar="character"),
          make_option(c("-s","--samples"), type="character", default=NULL, help="A comma-separated string of samples to plot individual samples"),
            make_option(c("-t","--trio"), action="store_true", default=FALSE, dest="trio", help="Boolean flag to denote if dataset is for trios [default = false]")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#
##
InputDir = opt$input
OutputDir <- paste0(InputDir,"plots")
datSet = opt$dataset
filter = opt$filter
samples = opt$samples
IndivDir <- paste0(InputDir,"indiv_dir")
trio = opt$trio
##
#
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input directory must be supplied", call.=FALSE)
}
#
if (is.null(opt$dataset)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (dataset name)", call.=FALSE)
}
#
#if (is.null(opt$filter)){
#  ifelse(!file.exists(file.path(InputDir,"filter_samples.txt")), file.create(file.path(InputDir,"filter_samples.txt")), FALSE)
#}
#
## DEBUG AREA ##
print(InputDir)
print(OutputDir)
print(datSet)
print(filter)
print(IndivDir)
print(trio)
#
stop("Execution halted. Debug mode active", call.=FALSE)
##
#
## Here we change to the working directory and make sure that the folder structure is correct and all of the files exist as input.
#
ifelse(!dir.exists(file.path(OutputDir)), dir.create(file.path(OutputDir)), FALSE)
ifelse(!dir.exists(file.path(IndivDir)), dir.create(file.path(IndivDir)), FALSE)
#
#datSet <- "" # Used if using something like Rstudio to directly run the script
TRIO=FALSE

samp.rm <- read.table(opt$filter,
                     sep = "\t")[,1] # If there are no samples to remove, simply leave file blank.

file1 <- paste(InputDir,datSet,".SNPs.AB_ref.all.txt", sep= "")
dat.ab_ref <- read.table(file1,
                     sep = "\t", col.names = c("SAMPLE", "AB", "COUNT"))
#samp <- unique(dat.ab_ref$SAMPLE) #deprecated line but if the rm_list.txt file is missing, this line can be substituted for not having samp.rm
samp <- samp[!(samp %in% samp.rm)]
dat.ab_ref <- dat.ab_ref[dat.ab_ref$SAMPLE %in% samp,]
ab_ref <- unique(dat.ab_ref$AB)
ab_ref.mat <- matrix(dat.ab_ref$COUNT,nrow = length(samp),ncol = length(ab_ref),
                 dimnames = list(samp,ab_ref), byrow = TRUE)

write.table(samp,
            paste(OutputDir,datSet,"_samps_all_used.txt", sep= ""),
            sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)

## This section now will import data frames for each recorded data point, allele balance (AB), genotype quality (GQ), and genotype (GT). ##

file1 <- paste(opt$input,datSet,".SNPs.AB_het.all.txt", sep= "")
dat.ab_het <- read.table(file1,
                  sep = "\t", col.names = c("SAMPLE", "AB", "COUNT"))
dat.ab_het <- dat.ab_het[dat.ab_het$SAMPLE %in% samp,]
ab_het <- unique(dat.ab_het$AB)
ab_het.mat <- matrix(dat.ab_het$COUNT,nrow = length(samp),ncol = length(ab_het),
                 dimnames = list(samp,ab_het), byrow = TRUE)


file1 <- paste(opt$input,datSet,".SNPs.AB_alt.all.txt", sep= "")
dat.ab_alt <- read.table(file1,
                         sep = "\t", col.names = c("SAMPLE", "AB", "COUNT"))
dat.ab_alt <- dat.ab_alt[dat.ab_alt$SAMPLE %in% samp,]
ab_alt <- unique(dat.ab_alt$AB)
ab_alt.mat <- matrix(dat.ab_alt$COUNT,nrow = length(samp),ncol = length(ab_alt),
                     dimnames = list(samp,ab_alt), byrow = TRUE)


file1 <- paste(opt$input,datSet,".SNPs.DP.all.txt", sep= "")
dat.dp <- read.table(file1,
                     sep = "\t", col.names = c("SAMPLE", "DP", "COUNT"))
dat.dp <- dat.dp[dat.dp$SAMPLE %in% samp,]
dp <- unique(dat.dp$DP)
dp.mat <- matrix(dat.dp$COUNT,nrow = length(samp),ncol = length(dp),
                 dimnames = list(samp,dp), byrow = TRUE)

file1 <- paste(opt$input,datSet,".SNPs.GQ_ref.all.txt", sep= "")
dat.gq_ref <- read.table(file1,
                         sep = "\t", col.names = c("SAMPLE", "GQ", "COUNT"))
dat.gq_ref <- dat.gq_ref[dat.gq_ref$SAMPLE %in% samp,]
gq_ref <- unique(dat.gq_ref$GQ)
gq_ref.mat <- matrix(dat.gq_ref$COUNT,nrow = length(samp),ncol = length(gq_ref),
                     dimnames = list(samp,gq_ref), byrow = TRUE)

file1 <- paste(opt$input,datSet,".SNPs.GQ_het.all.txt", sep= "")
dat.gq_het <- read.table(file1,
                         sep = "\t", col.names = c("SAMPLE", "GQ", "COUNT"))
dat.gq_het <- dat.gq_het[dat.gq_het$SAMPLE %in% samp,]
gq_het <- unique(dat.gq_het$GQ)
gq_het.mat <- matrix(dat.gq_het$COUNT,nrow = length(samp),ncol = length(gq_het),
                     dimnames = list(samp,gq_het), byrow = TRUE)

file1 <- paste(opt$input,datSet,".SNPs.GQ_alt.all.txt", sep= "")
dat.gq_alt <- read.table(file1,
                         sep = "\t", col.names = c("SAMPLE", "GQ", "COUNT"))
dat.gq_alt <- dat.gq_alt[dat.gq_alt$SAMPLE %in% samp,]
gq_alt <- unique(dat.gq_alt$GQ)
gq_alt.mat <- matrix(dat.gq_alt$COUNT,nrow = length(samp),ncol = length(gq_alt),
                     dimnames = list(samp,gq_alt), byrow = TRUE)

file1 <- paste(opt$input,datSet,".SNPs.GT.all.txt", sep= "")
dat.gt <- read.table(file1,
                     sep = "\t", col.names = c("SAMPLE", "GT", "COUNT"))
dat.gt <- dat.gt[dat.gt$SAMPLE %in% samp,]
gt <- unique(dat.gt$GT)
gt.mat <- matrix(dat.gt$COUNT,nrow = length(samp),ncol = length(gt),
                     dimnames = list(samp,gt), byrow = TRUE)


file1 <- paste(opt$input,datSet,".SNPs.GT.sum.all.txt", sep= "")
##
#
dat.gt_sum <- read.table(file1,
                         sep = "\t", col.names = c("SAMPLE", "GT", "COUNT", "REGION"))
dat.gt_sum <- dat.gt_sum[dat.gt_sum$SAMPLE %in% samp,]
gt <- unique(dat.gt_sum$GT)
region <- unique(dat.gt_sum$REGION)

dat.gt_sum$REGION <- factor(dat.gt_sum$REGION, levels = region)


gt_region.list <- list()
length(gt_region.list) <- length(samp)
names(gt_region.list) <- samp

for(i in names(gt_region.list)){
  dogSel <- i
  
  dog.df <- dat.gt_sum[dat.gt_sum$SAMPLE == dogSel,]
  dog.df.reg <- data.frame(ref = data.frame(dog.df[dog.df$GT == "0/0","COUNT"], row.names = dog.df[dog.df$GT == "0/1","REGION"])[region,],
                           het = data.frame(dog.df[dog.df$GT == "0/1","COUNT"], row.names = dog.df[dog.df$GT == "0/1","REGION"])[region,],
                           alt = data.frame(dog.df[dog.df$GT == "1/1","COUNT"], row.names = dog.df[dog.df$GT == "1/1","REGION"])[region,],
                           mis = data.frame(dog.df[dog.df$GT == "./.","COUNT"], row.names = dog.df[dog.df$GT == "./.","REGION"])[region,],
                           row.names = region)
  
  gt_region.list[[i]] <- dog.df.reg
  
}

# this is the GTs for every 5MB for every sample

gt_count_region.mat <- matrix(unlist(lapply(gt_region.list, rowSums, na.rm = TRUE)),
                              nrow = length(region), ncol = length(samp), dimnames = list(region,samp))

gt_Z_region.mat <- t(scale(t(gt_count_region.mat)))



# now to integrate into QC


# need to go through and define criteria for samples to be included in analysis
# criteria needs to be fairly strict.

# one way will be to go through and look at each metric as a histogram and pick outliers
# then choose plots that are generally agreed upon
# show which points are the outliers

# median coverage by missingness?




# samples were initially filtered according to missingness and coverage
# 10% missingness


dp.mean.all <- NULL
for(i in samp){
  dp.mean.all <- c(dp.mean.all, weighted.mean(0:99,w=dp.mat[i,]))
}


dp.cov10.all <- dp.cov5.all <- NULL
for(i in samp){
  dp.cov10.all <- c(dp.cov10.all, sum(dp.mat[i,10:ncol(dp.mat)])/sum(dp.mat[i,]))
  dp.cov5.all <- c(dp.cov5.all, sum(dp.mat[i,5:ncol(dp.mat)])/sum(dp.mat[i,]))
  
}


missingness.all <- gt.mat[,"./."]/rowSums(gt.mat, na.rm = TRUE) * 100


smoothScatter(x = dp.mean.all,
              y = log10(missingness.all))
plot(x = dp.mean.all,
     y = (missingness.all))
abline(v = 15, h = -3)

cols <- rep(1,length(missingness.all))
#cols[log10(missingness.all)>-3] <- 2
plot(dp.cov10.all,missingness.all)
plot(dp.mean.all[cols==1],dp.cov10.all[cols==1], col = scales::alpha(1,.5), pch = 16)
points(dp.mean.all[cols==2],dp.cov10.all[cols==2], col = scales::alpha(2,1), pch = 16)


## This section will isolate specific groups of samples and highlight them. ##
#
if(datSet == "CoolDataset"){
pdf(file = paste(OutputDir,datSet,".cov_10x_plots.coyote.pdf", sep= ""))
plot(dp.mean.all,dp.cov10.all, col = 1, pch = 16, cex = .5,
     xlab = "mean coverage", ylab = ">10X coverage", main = datSet)
points(dp.mean.all[grep("CLATUS",samp)],
       dp.cov10.all[grep("CLATUS",samp)], col = 2, pch = 16)
text(x = dp.mean.all[grep("CLATUS",samp)], y = dp.cov10.all[grep("CLATUS",samp)], 
     labels = samp[grep("CLATUS",samp)], pos = 4)
dev.off()
}


## Here plots samples with a mean coverage => 10X ##

pdf(file = paste(OutputDir,datSet,".cov_10x_plots.pdf", sep= ""),
    onefile = TRUE)
plot(dp.mean.all, dp.cov10.all, 
     xlab = "mean coverage", ylab = ">10X coverage", main = datSet,
     pch = 16, cex = .5)
plot(dp.mean.all, dp.cov10.all, 
     xlab = "mean coverage", ylab = ">10X coverage", main = datSet,
     pch = 16, cex = .5)
abline(v = 15, col = 2, lty = 2)
abline(h = .95, col = 2, lty = 2)
dev.off()



# next we analysed the per sample, per genotype distributions for allele balance (AB) and genotype quality (GQ)


# define cutoff here

plot(log10(colSums(ab_alt.mat[dp.mean.all >= 15,], na.rm = TRUE)))
plot((colSums(ab_alt.mat, na.rm = TRUE)))

plot((colSums(ab_het.mat[dp.mean.all >= 15,], na.rm = TRUE)))
plot((colSums(ab_ref.mat, na.rm = TRUE)))



samp.all <- samp
#samp <- samp[dp.mean.all >= 15] 

if(TRIO){
trio_samps <- read.table(opt$input,"trio_samples.txt", sep= "\t", header = TRUE)
samp <- samp[samp %in% trio_samps$Sample]
datSet <- "Trio"

}else{samp <- samp[dp.mean.all >= 15]}


ab_het.mat <- ab_het.mat[samp,]
ab_ref.mat <- ab_ref.mat[samp,]
ab_alt.mat <- ab_alt.mat[samp,]

gq_het.mat <- gq_het.mat[samp,]
gq_ref.mat <- gq_ref.mat[samp,]
gq_alt.mat <- gq_alt.mat[samp,]

dp.mat <- dp.mat[samp,]

gt.mat <- gt.mat[samp,]

names(dp.mean.all) <- samp.all
names(dp.cov10.all) <- samp.all
names(dp.cov5.all) <- samp.all

dp.mean <- dp.mean.all[samp]
dp.cov10 <- dp.cov10.all[samp]
dp.cov5 <- dp.cov5.all[samp]


gt_region.list <- gt_region.list[samp]
gt_count_region.mat <- gt_count_region.mat[,samp]
gt_Z_region.mat <- gt_Z_region.mat[,samp]

  
outlier.cut <- abs(gt_Z_region.mat) > 10
outlier_count_regions.df <- data.frame(SAMPLE = colnames(gt_Z_region.mat)[col(gt_Z_region.mat)[outlier.cut]],
                                       REGION = rownames(gt_Z_region.mat)[row(gt_Z_region.mat)[outlier.cut]],
                                       Z = gt_Z_region.mat[outlier.cut],
                                       GTs = gt_count_region.mat[outlier.cut])

outlier_count_regions.df <- outlier_count_regions.df[outlier_count_regions.df$GTs == 0,]  

# counted genotype annotations
# does any sample have excess non-counted genotype annotations?

sd.lim <- c(mean(rowSums(gt.mat)) - (5*sd(rowSums(gt.mat))),
            mean(rowSums(gt.mat)) + (5*sd(rowSums(gt.mat))))
x.lim <- range(c(range(rowSums(gt.mat)),sd.lim))

pdf(file = paste(OutputDir,datSet,".counted_annotations.pdf", sep= ""),
    onefile = TRUE)

# depleted annotations

hist(rowSums(gt.mat), breaks = 100, xlim = x.lim,
     border = NA, col = 8,
     xlab = "Counted GT annotations",
     main = datSet)
abline(v = mean(rowSums(gt.mat)), col = 2, lwd = 1)
abline(v = mean(rowSums(gt.mat)) - (5*sd(rowSums(gt.mat))), col = 2,lty = 2, lwd = 1)
abline(v = mean(rowSums(gt.mat)) + (5*sd(rowSums(gt.mat))), col = 2,lty = 2, lwd = 1)

plot(rowSums(gt_Z_region.mat < -10), type = "h",
     xlab = "Genomic region",
     ylab = "Samples with excess non-counted GT annotations",
     las = 1, 
     main = datSet)

plot(rowSums(gt_count_region.mat == 0), type = "h",
     xlab = "Genomic region",
     ylab = "Samples with no counted GT annotations",
     las = 1, 
     main = datSet)


dev.off()

write.table(outlier_count_regions.df,
            file = paste(OutputDir,datSet,".outlier_GTcount_regions.txt", sep= ""),
            sep = "\t",quote = FALSE,row.names = FALSE, col.names = TRUE)

  
if(TRIO){
pdf(file = paste(OutputDir,datSet,".cov_10x_plots.pdf", sep= ""))
plot(dp.mean, dp.cov10, 
     xlab = "mean coverage", ylab = ">10X coverage", main = datSet,
     pch = 16, cex = .5)
dev.off()
}

yVals <- rollmean(colSums(ab_het.mat, na.rm = TRUE),k = 10,na.pad = TRUE)
yVals[is.na(yVals)] <- 0
yVals <- yVals/sum(yVals) * 100


pdf(file = paste(OutputDir,datSet,".filter_plots.pdf", sep= ""),
    onefile = TRUE)


plot(yVals,
     type = "l", 
     xlab = "Hetrozygous allele balance (%)",
     ylab = "Genotypes (%)")

segments(x0 = seq(10,30,5),
         x1 = seq(10,30,5),
         y0 = 0,
         y1 = yVals[seq(10,30,5)],
         lty = 2)
segments(x0 = seq(70,90,5),
         x1 = seq(70,90,5),
         y0 = 0,
         y1 = yVals[seq(70,90,5)],
         lty = 2)

segments(x0 = 0,
         x1 = seq(10,30,5),
         y0 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE),
         y1 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE))
segments(x0 = seq(10,30,5),
         x1 = seq(10,30,5),
         y0 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE) - .01*max(yVals,na.rm = TRUE),
         y1 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE))
segments(x0 = 0,
         x1 = 0,
         y0 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE) - .01*max(yVals,na.rm = TRUE),
         y1 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE))

segments(x0 = 100,
         x1 = rev(seq(70,90,5)),
         y0 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE),
         y1 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE))
segments(x0 = rev(seq(70,90,5)),
         x1 = rev(seq(70,90,5)),
         y0 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE) - .01*max(yVals,na.rm = TRUE),
         y1 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE))
segments(x0 = 100,
         x1 = 100,
         y0 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE) - .01*max(yVals,na.rm = TRUE),
         y1 = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE))
text(x = seq(10,30,5)/2,
     y = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE),
     labels = paste(seq(10,30,5),"%", sep = ""), c(0.5,-0.1))
text(x = rev(seq(70,90,5)) + seq(10,30,5)/2,
     y = yVals[seq(10,30,5)] + seq(.02,.06,.01) * max(yVals,na.rm = TRUE),
     labels = paste(seq(10,30,5),"%", sep = ""), c(0.5,-0.1))



# how many uncertain hets are there?
het_tail_30 <- het_tail_25 <- het_tail_20 <- het_tail_15 <- het_tail_10 <- NULL
for(i in samp){
het_tail_30 <- c(het_tail_30,sum(ab_het.mat[i,c(1:30,71:100)])/sum(ab_het.mat[i,]))
het_tail_25 <- c(het_tail_25,sum(ab_het.mat[i,c(1:25,76:100)])/sum(ab_het.mat[i,]))
het_tail_20 <- c(het_tail_20,sum(ab_het.mat[i,c(1:20,81:100)])/sum(ab_het.mat[i,]))
het_tail_15 <- c(het_tail_15,sum(ab_het.mat[i,c(1:15,86:100)])/sum(ab_het.mat[i,]))
het_tail_10 <- c(het_tail_10,sum(ab_het.mat[i,c(1:10,91:100)])/sum(ab_het.mat[i,]))
}


# how many samples have less than x
plot(sort(het_tail_30) * 100, (1:length(samp))/length(samp) * 100, type = "l",
     ylab = "Samples (percentile)", xlab = "Failed het genotypes (%)",
     xlim = c(0,max(het_tail_30)*100))
lines(sort(het_tail_25) * 100, (1:length(samp))/length(samp) * 100, xlim = c(0,1), type = "l")
lines(sort(het_tail_20) * 100, (1:length(samp))/length(samp) * 100, xlim = c(0,1), type = "l")
lines(sort(het_tail_15) * 100, (1:length(samp))/length(samp) * 100, xlim = c(0,1), type = "l")
lines(sort(het_tail_10) * 100, (1:length(samp))/length(samp) * 100, xlim = c(0,1), type = "l")

points(x = tail(sort(het_tail_30)[sort(het_tail_30)<=.05], n = 1) * 100,
       y = tail((1:length(samp))[sort(het_tail_30)<=.05], n = 1)/length(samp) * 100,
       pch = 16, col = 2)

points(x = tail(sort(het_tail_25)[sort(het_tail_25)<=.05], n = 1) * 100,
       y = tail((1:length(samp))[sort(het_tail_25)<=.05], n = 1)/length(samp) * 100,
       pch = 16, col = 2)

points(x = tail(sort(het_tail_20)[sort(het_tail_20)<=.05], n = 1) * 100,
       y = tail((1:length(samp))[sort(het_tail_20)<=.05], n = 1)/length(samp) * 100,
       pch = 16, col = 2)

points(x = tail(sort(het_tail_15)[sort(het_tail_15)<=.05], n = 1) * 100,
       y = tail((1:length(samp))[sort(het_tail_15)<=.05], n = 1)/length(samp) * 100,
       pch = 16, col = 2)

points(x = tail(sort(het_tail_10)[sort(het_tail_10)<=.05], n = 1) * 100,
       y = tail((1:length(samp))[sort(het_tail_10)<=.05], n = 1)/length(samp) * 100,
       pch = 16, col = 2)



text(x = tail(sort(het_tail_30)[sort(het_tail_30)<=.05], n = 1) * 100,
     y = tail((1:length(samp))[sort(het_tail_30)<=.05], n = 1)/length(samp) * 100,
     labels = " 30 %",adj = c(0,1))

text(x = tail(sort(het_tail_25)[sort(het_tail_25)<=.05], n = 1) * 100,
     y = tail((1:length(samp))[sort(het_tail_25)<=.05], n = 1)/length(samp) * 100,
     labels = " 25 %",adj = c(0,1))

text(x = tail(sort(het_tail_20)[sort(het_tail_20)<=.05], n = 1) * 100,
     y = tail((1:length(samp))[sort(het_tail_20)<=.05], n = 1)/length(samp) * 100,
     labels = " 20 %",adj = c(0,1))

text(x = tail(sort(het_tail_15)[sort(het_tail_15)<=.05], n = 1) * 100,
     y = tail((1:length(samp))[sort(het_tail_15)<=.05], n = 1)/length(samp) * 100,
     labels = " 15 %",adj = c(1,0))

text(x = tail(sort(het_tail_10)[sort(het_tail_10)<=.05], n = 1) * 100,
     y = tail((1:length(samp))[sort(het_tail_10)<=.05], n = 1)/length(samp) * 100,
     labels = " 10 %",adj = c(0,0))


# other version of this plot

plot(sort(het_tail_30) * 100, (1:length(samp))/length(samp) * 100, type = "l",
     ylab = "Samples (percentile)", xlab = "Failed het genotypes (%)",
     xlim = c(0,max(het_tail_30)*100), lty = 2)
lines(sort(het_tail_25) * 100, (1:length(samp))/length(samp) * 100, xlim = c(0,1), type = "l", lty = 3)
lines(sort(het_tail_20) * 100, (1:length(samp))/length(samp) * 100, xlim = c(0,1), type = "l", lty = 4)
lines(sort(het_tail_15) * 100, (1:length(samp))/length(samp) * 100, xlim = c(0,1), type = "l", lty = 5)
lines(sort(het_tail_10) * 100, (1:length(samp))/length(samp) * 100, xlim = c(0,1), type = "l", lty = 6)

abline(v = 5, lty = 1, col = 2)

legend("right",
       # legend = c("30%, 70%",
       #            "25%, 75%",
       #            "20%, 80%",
       #            "15%, 85%",
       #            "10%, 90%"),
       legend = c("30%",
                  "25%",
                  "20%",
                  "15%",
                  "10%"),
       lty = 2:6,
       title = "AB cutoffs", 
       bty = "n")



hist(het_tail_20 * 100, breaks = seq(0,100,.5), 
     xlab = "Outlier variants (%)",
     ylab = "Samples (n)",
     main = paste(datSet,"AB balance (0/1)"),
     xlim = c(0,40), border = NA, col =8)
abline(v = 5, col = 2)




plot((colSums(ab_alt.mat, na.rm = TRUE)), type = "h",
     main = paste(datSet,"AB 1/1"),
     ylab = "Variants (n)", xlab = "Ref read (%)");abline(v = 5, col = 2)
plot((colSums(ab_ref.mat, na.rm = TRUE)),type = "h",
      main = paste(datSet,"AB 0/0"),
      ylab = "Variants (n)", xlab = "Ref read (%)");abline(v = 95, col = 2)

# we could put points in on the number of samples
# keep the cutoff at 5%, what boundaries can we have?
# plot multiple lines for each cut off range

ref_tail<- NULL
for(i in samp){
  ref_tail <- c(ref_tail,sum(ab_ref.mat[i,c(1:94)])/sum(ab_ref.mat[i,]))
}
hist(ref_tail * 100, breaks = seq(0,100,1), 
     xlab = "Outlier variants (%)",
     ylab = "Samples",
     main = "AB balance (0/0)",
     border = NA, col =8)
abline(v = 5, col = 2)


alt_tail<- NULL
for(i in samp){
  alt_tail <- c(alt_tail,1 - (sum(ab_alt.mat[i,c(1:5)])/sum(ab_alt.mat[i,])))
}
hist(alt_tail * 100, breaks = seq(0,100,1), 
     xlab = "Outlier variants (%)",
     ylab = "Samples",
     main = "AB balance (1/1)",
     border = NA, col = 8)
abline(v = 5, col = 2)


v.diag <- venn::venn(list(samp[het_tail_20 > .05],
                          samp[ref_tail > .05],
                          samp[alt_tail > .05]),
                     snames = c("0/1","0/0","1/1"),
                     main = paste(datSet,"filtered AB"),
                     par = FALSE, box = FALSE)


# so now we have the GQ stuff to filter on


 plot((colSums(gq_alt.mat)), type = "h")
 plot((colSums(gq_het.mat)), type = "h")
 plot((colSums(gq_ref.mat)), type = "h")
 
# 
 plot(cumsum(colSums(gq_alt.mat))/sum(gq_alt.mat), type = "h")
 abline(h = .05)
# 
# # at which point to we capture the bottom 5% of sites
# 
 plot(cumsum(colSums(gq_ref.mat))/sum(gq_ref.mat), type = "h")
 abline(h = .05)
# 
 plot(cumsum(colSums(gq_het.mat))/sum(gq_het.mat), type = "h")
 abline(h = .05)


# 30 30 80

# samples enriched for the bottom 5% of sites

gq_alt.tail <- rowSums(gq_alt.mat[,1:20])/rowSums(gq_alt.mat)
gq_ref.tail <- rowSums(gq_ref.mat[,1:20])/rowSums(gq_ref.mat)
gq_het.tail <- rowSums(gq_het.mat[,1:20])/rowSums(gq_het.mat)

plot(data.frame(gq_alt.tail,gq_het.tail,gq_ref.tail))
abline(h = .05)



plot((0:99)-.3, 
     colSums(gq_ref.mat)/sum(gq_ref.mat), 
     type = "h", ylim = c(0,.2), 
     xlab = "GQ", ylab = "Proportion",
     lwd = 2, main = datSet)
points(0:99,
       colSums(gq_het.mat)/sum(gq_het.mat), 
       type = "h", col = 2, lwd = 2)
points((0:99)+.3, 
       colSums(gq_alt.mat)/sum(gq_alt.mat), 
       type = "h", col = 4, lwd = 2)
legend("topleft", legend = c("0/0", "0/1", "1/1"),
       col = c(1,2,4), lty = 1, lwd = 2, 
       title = "Genotype", bty = "n")


# samples and failed GT
hist(gq_het.tail*100, breaks = 100, 
     xlab = "Outlier variants (%)",
     ylab = "Samples (n)",
     main = paste(datSet,"GQ (0/1)"),
     border = NA, col = 8); abline(v = 5, col = 2)
hist(gq_alt.tail*100, breaks = 100, 
     xlab = "Outlier variants (%)",
     ylab = "Samples (n)",
     main = paste(datSet,"GQ (1/1)"),
     border = NA, col = 8); abline(v = 5, col = 2)
hist(gq_ref.tail*100, breaks = 100, 
     xlab = "Outlier variants (%)",
     ylab = "Samples (n)",
     main = paste(datSet,"GQ (0/0)"),
     border = NA, col = 8); abline(v = 5, col = 2)


v.diag <- venn::venn(list(samp[gq_het.tail > .05],
                          samp[gq_ref.tail > .05],
                          samp[gq_alt.tail > .05]),
                     snames = c("0/1","0/0","1/1"),
                     main = paste(datSet,"filtered GQ"),
                     par = FALSE, box = FALSE)

dev.off()



# confidence levels
# Phred 20 , 10^(-20/10) = 0.01, 1% chance genotype is incorrect?
# or is it 100 times more confidence?

# put together a big table of all the criteria
filter_rm <- data.frame(ref.ab = ref_tail > 0.05,
                        het.ab = het_tail_20 > 0.05,
                        alt.ab = alt_tail > 0.05,
                        ref.gq = gq_ref.tail > 0.05,
                        het.gq = gq_het.tail > 0.05,
                        alt.gq = gq_alt.tail > 0.05,
                        cov.10 = dp.cov10 < 0.95,
                        count.gt = rowSums(gt.mat) < sd.lim[1],
                        region.gt = samp %in% unique(outlier_count_regions.df$SAMPLE))
rownames(filter_rm) <- samp

filter_rm.sort <- filter_rm[order(dp.mean),]


pdf(file = paste(OutputDir,datSet,".all_filters.pdf", sep= ""))

layout(1:2)
par(mar = c(0,5,4,1), oma = c(5,5,5,1))

plot((1:length(samp))-.5,
     sort(dp.mean), 
     type = "h",
     xaxs = "i", yaxs = "i",
     xlim = c(0,length(samp)),
     ylim = c(0,max(dp.mean)+5),
     ylab = "Coverage", xaxt = "n", 
     bty = "n", xlab = "", las = 2,
     main = datSet)

par(mar = c(2,5,1,1))

plot(0,0, type = "n",
     xlim = c(0,length(samp)),
     ylim = c(-1.5,ncol(filter_rm.sort)),
     xaxs = "i", yaxs = "i", axes = FALSE,
     xlab = "", ylab = "")

for(i in 1:ncol(filter_rm.sort)){
  select_filt <- (1:nrow(filter_rm.sort))[filter_rm.sort[,i]]
  if(length(select_filt) > 0){
  rect(xleft = select_filt - 1,
       xright = select_filt,
       ybottom = i - 1,
       ytop = i, 
       col = "black", 
       density = -1, border = FALSE)
  }
}

sum.cols <- RColorBrewer::brewer.pal(7,"GnBu")

failedSamp <- rowSums(filter_rm.sort)
failedSamp.no <- which(failedSamp > 0)
failedSamp.col <- sum.cols[failedSamp[failedSamp>0]]

rect(xleft = failedSamp.no - 1,
     xright = failedSamp.no,
     ytop = -.5,
     ybottom = -1.5,
     col = failedSamp.col, density = -1,
     border = FALSE)

mtext(side = 2, at = 1:ncol(filter_rm.sort) - .5,
      text = colnames(filter_rm.sort), las = 2, line = .5)
mtext(side = 2, at = -.75,
      text = "Total", las = 2, line = .5)

dev.off()



pdf(file = paste(OutputDir,datSet,".upset.pdf", sep= ""))
filter_upset <- apply(filter_rm.sort, MARGIN = 2,FUN = as.integer)
colnames(filter_upset) <- colnames(filter_rm.sort)
upset(data.frame(filter_upset), nsets = 9,
      point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Filter Intersections", sets.x.label = "Filtered Samples",
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), order.by = "freq")
dev.off()

write.table(x = data.frame(SAMP = rownames(filter_rm), 
                       filter_rm, 
                       KEEP = !(rowSums(filter_rm) > 0)),
            sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE,
            file = paste(OutputDir,datSet,".filters.txt", sep= ""))

seq_stats <- data.frame(ref.ab = ref_tail,
                        het.ab = het_tail_20,
                        alt.ab = alt_tail,
                        ref.gq = gq_ref.tail,
                        het.gq = gq_het.tail,
                        alt.gq = gq_alt.tail,
                        cov.10 = dp.cov10,
                        cov = dp.mean[samp],
                        KEEP = rowSums(filter_rm) == 0,
                        gt.mat[samp,],
                        gt.count = rowSums(gt.mat[samp,]),
                        fail.gt.count.regions = table(factor(outlier_count_regions.df$SAMPLE, levels = samp))[samp])
rownames(seq_stats) <- samp
colnames(seq_stats)[(ncol(seq_stats)-6):ncol(seq_stats)] <- c("0/0","0/1","1/1","./.","GT Count","Fail Sample ID","Fail Regions count")

write.table(x = data.frame(SAMP = rownames(seq_stats), 
                           seq_stats),
            sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE,
            file = paste(OutputDir,datSet,".seq_stats.txt", sep= ""))


colSums(filter_rm)
filter_rm[filter_rm$region.gt,]



sum(rowSums(filter_rm) != 0)

sum(rowSums(filter_rm) == 0)


dim(gt.mat)


cols <- rep(2, length(samp.all))
cols[samp.all %in% rownames(filter_rm)[rowSums(filter_rm) < 1]] <- 1
plot(gt.mat[,2]/rowSums(gt.mat[,c(1,3)]), col = cols)



# plots for individual samples compared to whole distribution, will be edited in, tested, and fixed later.

samples.defined = function(samples)!is.null(samples)
#
if(samples.defined){
  indiv.samples <- strsplit(samples, ",")[[1]]
  #indiv.samples <- trimws(indiv.samples)


  i = indiv.samples[2]
  for(i in indiv.samples){
    pdf(file = paste(IndivDir,datSet,".",i,".stats.pdf", sep= ""),
        onefile = TRUE)
    
  plot(dp.mean.all,dp.cov10.all, col = 1, pch = 16, cex = .5,
       xlab = "mean coverage", ylab = ">10X coverage", main = datSet)
  points(dp.mean.all[i],
         dp.cov10.all[i], col = 2, pch = 16)
  text(x = dp.mean.all[i], y = dp.cov10.all[i], 
       labels = i, pos = 1, col= 2)
  
  
  plot((0:99)-.25, 
       colSums(gq_ref.mat)/sum(gq_ref.mat), 
       type = "h", ylim = c(0,.2), 
       xlab = "GQ Ref", ylab = "Proportion",
       lwd = 2, main = paste(datSet,"GQ (0/0)"))
  points((0:99)+.25,
         gq_ref.mat[i,]/sum(gq_ref.mat[i,]),
         type = "h", col = 2, lwd = 2)
  legend("topright", lwd = c(2,2), col = c(1,2), legend = c("All Dogs", i))
  
  plot((0:99)-.25, 
       colSums(gq_het.mat)/sum(gq_het.mat), 
       type = "h", ylim = c(0,.2), 
       xlab = "GQ Het", ylab = "Proportion",
       lwd = 2, main = paste(datSet,"GQ (0/1)"))
  points((0:99)+.25,
         gq_het.mat[i,]/sum(gq_het.mat[i,]),
         type = "h", col = 2, lwd = 2)
  legend("topright", lwd = c(2,2), col = c(1,2), legend = c("All Dogs", i))
  
  plot((0:99)-.25, 
       colSums(gq_alt.mat)/sum(gq_alt.mat), 
       type = "h", ylim = c(0,.2), 
       xlab = "GQ Alt", ylab = "Proportion",
       lwd = 2, main = paste(datSet,"GQ (1/1)"))
  points((0:99)+.25,
         gq_alt.mat[i,]/sum(gq_alt.mat[i,]),
         type = "h", col = 2, lwd = 2)
  #legend("topright", lwd = c(2,2), col = c(1,2), legend = c("All Dogs", i))
  
  
  
  hist(gq_het.tail*100, breaks = 100, 
       xlab = "Outlier variants (%)",
       ylab = "Samples (n)",
       main = paste(datSet,"GQ (0/1)")); abline(v = 5, col = 2, lty = 2, lwd = 2) ; abline(v = (gq_het.tail*100)[i], col = 4, lwd = 2)
  legend("topright", lwd = c(2,2), lty = c(2,1), col = c(2,4), legend = c("Cutoff", i))
  hist(gq_alt.tail*100, breaks = 100, 
       xlab = "Outlier variants (%)",
       ylab = "Samples (n)",
       main = paste(datSet,"GQ (1/1)")); abline(v = 5, col = 2, lty = 2, lwd = 2) ; abline(v = (gq_alt.tail*100)[i], col = 4, lwd = 2)
  legend("topright", lwd = c(2,2), lty = c(2,1), col = c(2,4), legend = c("Cutoff", i))
  hist(gq_ref.tail*100, breaks = 100, 
       xlab = "Outlier variants (%)",
       ylab = "Samples (n)",
       main = paste(datSet,"GQ (0/0)")); abline(v = 5, col = 2, lty = 2, lwd = 2) ; abline(v = (gq_ref.tail*100)[i], col = 4, lwd = 2)
  legend("topright", lwd = c(2,2), lty = c(2,1), col = c(2,4), legend = c("Cutoff", i))
  
  
  
  plot((0:99)-.25,
       (colSums(ab_het.mat, na.rm = TRUE))/sum(ab_het.mat,na.rm = TRUE) * 100, type = "h",
       main = paste(datSet,"AB 1/0"),
       ylab = "Variants (%)", xlab = "Ref read (%)", lwd = 2)
  points((0:99)+.25,
         ab_het.mat[i,]/sum(ab_het.mat[i,]) * 100, type = "h", col = 2, lwd = 2)
  legend("topright", lwd = c(2,2), col = c(1,2), legend = c("All Dogs", i))
  
  plot((0:99)-.25,
       (colSums(ab_alt.mat, na.rm = TRUE))/sum(ab_alt.mat,na.rm = TRUE) * 100, type = "h",
       main = paste(datSet,"AB 1/1"),
       ylab = "Variants (%)", xlab = "Ref read (%)", lwd = 2)
  points((0:99)+.25,
         ab_alt.mat[i,]/sum(ab_alt.mat[i,]) * 100, type = "h", col = 2, lwd = 2)
  legend("topright", lwd = c(2,2), col = c(1,2), legend = c("All Dogs", i))
  
  plot((0:99)-.25,
       (colSums(ab_ref.mat, na.rm = TRUE))/sum(ab_ref.mat,na.rm = TRUE) * 100, type = "h",
       main = paste(datSet,"AB 0/0"),
       ylab = "Variants (%)", xlab = "Ref read (%)", lwd = 2)
  points((0:99)+.25,
         ab_ref.mat[i,]/sum(ab_ref.mat[i,]) * 100, type = "h", col = 2, lwd = 2)
  legend("topright", lwd = c(2,2), col = c(1,2), legend = c("All Dogs", i))
  
  
  
  plot(c(0,sort(het_tail_30) * 100,100), c(0,(1:length(samp))/length(samp) * 100,100), type = "l",
       ylab = "Samples (percentile)", xlab = "Failed het genotypes (%)",
       main = "Allele balance cutoff distributions",
       xlim = c(0,max(het_tail_30)*100))
  lines(c(0,sort(het_tail_25) * 100,100), c(0,(1:length(samp))/length(samp) * 100,100), xlim = c(0,1), type = "l")
  lines(c(0,sort(het_tail_20) * 100,100), c(0,(1:length(samp))/length(samp) * 100,100), xlim = c(0,1), type = "l")
  lines(c(0,sort(het_tail_15) * 100,100), c(0,(1:length(samp))/length(samp) * 100,100), xlim = c(0,1), type = "l")
  lines(c(0,sort(het_tail_10) * 100,100), c(0,(1:length(samp))/length(samp) * 100,100), xlim = c(0,1), type = "l")
  
  abline(v = 5, lty = 2, lwd = 2, col = 2)
  
  points(x = het_tail_30[samp == i] * 100,
         y = which(samp[order(het_tail_30)] == i)/length(samp) * 100,
         pch = 16, col = 4)
  points(x = het_tail_25[samp == i] * 100,
         y = which(samp[order(het_tail_25)] == i)/length(samp) * 100,
         pch = 16, col = 4)
  points(x = het_tail_20[samp == i] * 100,
         y = which(samp[order(het_tail_20)] == i)/length(samp) * 100,
         pch = 16, col = 4)
  points(x = het_tail_15[samp == i] * 100,
         y = which(samp[order(het_tail_15)] == i)/length(samp) * 100,
         pch = 16, col = 4)
  points(x = het_tail_10[samp == i] * 100,
         y = which(samp[order(het_tail_10)] == i)/length(samp) * 100,
         pch = 16, col = 4)
  
  
  text(x = sort(het_tail_30)[(length(samp)*.7)] * 100,
       y = 70,
       labels = " 30%",adj = c(0,1))
  
  text(x = sort(het_tail_25)[(length(samp)*.6)] * 100,
       y = 60,
       labels = " 25%",adj = c(0,1))
  
  text(x = sort(het_tail_20)[(length(samp)*.5)] * 100,
       y = 50,
       labels = " 20%",adj = c(0,1))
  
  text(x = sort(het_tail_15)[(length(samp)*.4)] * 100,
       y = 40,
       labels = " 15%",adj = c(0,1))
  
  text(x = sort(het_tail_10)[(length(samp)*.3)] * 100,
       y = 30,
       labels = " 10%",adj = c(0,1))
  
  legend("bottomright", col = 4, pch = 16, legend = i)
  
  
  hist(het_tail_20 * 100, breaks = 100, 
       xlab = "Outlier variants (%)",
       ylab = "Samples (n)",
       main = paste(datSet,"AB (0/1)"))
  abline(v = 5, col = 2, lty = 2, lwd = 2) ; 
  abline(v = (het_tail_20*100)[which(samp == i)], col = 4, lwd = 2)
  legend("topright", lwd = c(2,2), lty = c(2,1), col = c(2,4), legend = c("Cutoff", i))
  
  hist(ref_tail * 100, breaks = 100, 
       xlab = "Outlier variants (%)",
       ylab = "Samples",
       main = paste(datSet,"AB (0/0)"))
  abline(v = 5, col = 2, lty = 2, lwd = 2) ; 
  abline(v = (ref_tail*100)[which(samp == i)], col = 4, lwd = 2)
  legend("topright", lwd = c(2,2), lty = c(2,1), col = c(2,4), legend = c("Cutoff", i))
  
  hist(alt_tail * 100, breaks = 100, 
       xlab = "Outlier variants (%)",
       ylab = "Samples",
       main = paste(datSet,"AB (1/1)"))
  abline(v = 5, col = 2, lty = 2, lwd = 2) ; 
  abline(v = (alt_tail*100)[which(samp == i)], col = 4, lwd = 2)
  legend("topright", lwd = c(2,2), lty = c(2,1), col = c(2,4), legend = c("Cutoff", i))
  
  
  hist(gt.mat[,"./."]/rowSums(gt.mat), breaks = 100,
       xlab = "./. calls (%)", ylab = "Samples (n)",
       main = paste(datSet,"Calls (./.)"))
  abline(v = gt.mat[i,"./."]/rowSums(gt.mat)[i],  col = 4, lwd = 2)
  legend("topright", lwd = c(2), lty = c(1), col = c(4), legend = c(i))
  
  
  hist(gt.mat[,"0/0"]/rowSums(gt.mat[,c("0/0","0/1","1/1")]) * 100, breaks = 100,
       xlab = "0/0 calls (%)", ylab = "Samples (n)",
       main = paste(datSet,"Calls (0/0)"))
  abline(v = gt.mat[i,"0/0"]/rowSums(gt.mat[,c("0/0","0/1","1/1")])[i]  * 100, col = 4, lwd = 2)
  legend("topright", lwd = c(2), lty = c(1), col = c(4), legend = c(i))
  
  
  hist(gt.mat[,"0/1"]/rowSums(gt.mat[,c("0/0","0/1","1/1")]) * 100, breaks = 100,
       xlab = "0/1 calls (%)", ylab = "Samples (n)",
       main = paste(datSet,"Calls (0/1)"))
  abline(v = gt.mat[i,"0/1"]/rowSums(gt.mat[,c("0/0","0/1","1/1")])[i]  * 100, col = 4, lwd = 2)
  legend("topright", lwd = c(2), lty = c(1), col = c(4), legend = c(i))
  
  
  hist(gt.mat[,"1/1"]/rowSums(gt.mat[,c("0/0","0/1","1/1")]) * 100, breaks = 100,
       xlab = "1/1 calls (%)", ylab = "Samples (n)",
       main = paste(datSet,"Calls (1/1)"))
  abline(v = gt.mat[i,"1/1"]/rowSums(gt.mat[,c("0/0","0/1","1/1")])[i]  * 100, col = 4, lwd = 2)
  legend("topright", lwd = c(2), lty = c(1), col = c(4), legend = c(i))
  
  
  hist(gt.mat[,"0/1"]/gt.mat[,"1/1"], breaks = 100,
       xlab = "het/alt ratio", ylab = "Samples (n)",
       main = paste(datSet,"Het/alt ratio"))
  abline(v = gt.mat[i,"0/1"]/gt.mat[i,"1/1"],  col = 4, lwd = 2)
  legend("topright", lwd = c(2), lty = c(1), col = c(4), legend = c(i))
  
  dev.off()
  
  }
}












filter_rm

dp.mean


# plot the distributions for samples with different coverages 

pdf(file = paste(OutputDir,datSet,".GQ_cov.pdf", sep= ""),
    width = 8, height = 10)
layout(matrix(1:6,ncol = 1))
par(mar = c(0,4,0,1),oma = c(5,5,5,1))
plot(colSums(gq_ref.mat[dp.mean > 15 & dp.mean <= 20 & rowSums(filter_rm) == 0,]), ylab = "", type = "h", lwd = 2, xaxt = "n", las = 2); legend("topleft", legend = c("cov = 20-25", paste(sum(dp.mean > 15 & dp.mean <= 20 & rowSums(filter_rm) == 0)), "Samples"), bty = "n")
plot(colSums(gq_ref.mat[dp.mean > 20 & dp.mean <= 25 & rowSums(filter_rm) == 0,]), ylab = "", type = "h", lwd = 2, xaxt = "n", las = 2); legend("topleft", legend = c("cov = 20-25", paste(sum(dp.mean > 20 & dp.mean <= 25 & rowSums(filter_rm) == 0)), "Samples"), bty = "n")
plot(colSums(gq_ref.mat[dp.mean > 25 & dp.mean <= 30 & rowSums(filter_rm) == 0,]), ylab = "", type = "h", lwd = 2, xaxt = "n", las = 2); legend("topleft", legend = c("cov = 25-30", paste(sum(dp.mean > 25 & dp.mean <= 30 & rowSums(filter_rm) == 0)), "Samples"), bty = "n")
plot(colSums(gq_ref.mat[dp.mean > 30 & dp.mean <= 35 & rowSums(filter_rm) == 0,]), ylab = "", type = "h", lwd = 2, xaxt = "n", las = 2); legend("topleft", legend = c("cov = 30-35", paste(sum(dp.mean > 30 & dp.mean <= 35 & rowSums(filter_rm) == 0)), "Samples"), bty = "n")
plot(colSums(gq_ref.mat[dp.mean > 35 & dp.mean <= 40 & rowSums(filter_rm) == 0,]), ylab = "", type = "h", lwd = 2, xaxt = "n", las = 2); legend("topleft", legend = c("cov = 35-40", paste(sum(dp.mean > 35 & dp.mean <= 40 & rowSums(filter_rm) == 0)), "Samples"), bty = "n")
plot(colSums(gq_ref.mat[dp.mean > 40 & dp.mean <= 45 & rowSums(filter_rm) == 0,]), ylab = "", type = "h", lwd = 2, xaxt = "n", las = 2); legend("topleft", legend = c("cov = 40-45", paste(sum(dp.mean > 40 & dp.mean <= 45 & rowSums(filter_rm) == 0)), "Samples"), bty = "n")
#plot(colSums(gq_ref.mat[dp.mean > 45 & dp.mean <= 50 & rowSums(filter_rm) == 0,]), ylab = "", type = "h", lwd = 2, xaxt = "n", las = 2); legend("topleft", legend = c("cov = 45-50", paste(sum(dp.mean > 45 & dp.mean <= 50 & rowSums(filter_rm) == 0)), "Samples"), bty = "n")
axis(1,at = seq(0,100,10))
mtext("GQ", side = 1,outer = TRUE, line = 3)
mtext("Variants (0/0)", side = 2,outer = TRUE, line = 1.7)
title(main = "0/0 GQ distribution by WGS coverage", outer = TRUE)
dev.off()
# maybe we should lower the het gq criteria?

## Trio section to fix
#if(TRIO){
#trio_samps <- read.table("~/Desktop/oneDrive/project_data/joint_coding_regions/Sample_info/trio_samples_2216g.txt", sep= "\t", header = TRUE)


#cols <- rep(1, length(samp))
#cols[samp %in% trio_samps$Sample] <- 2

#plot(dp.mean,dp.cov10, col = 1, pch = 16, cex = .5)
#points(dp.mean[cols == 2], dp.cov10[cols == 2],
#       col = 2, pch = 16)

# maybe if I show the fail plot 
# then show all of the dogs below it as an extra row



#pdf(file = "~/Desktop/oneDrive/project_data/joint_coding_regions/trios/trio_samples_2216g_failed_filters_only.pdf",
#    height = 7, width = 8)
#layout(1:3, height = c(.5,1,.5))
#par(mar = c(0,5,4,1), oma = c(8,5,5,1))

#plot((1:length(samp))-.5,
#     sort(dp.mean), 
#     type = "h",
#     xaxs = "i", yaxs = "i",
#     xlim = c(0,length(samp)),
#     ylim = c(0,max(dp.mean)+5),
#     ylab = "Coverage", xaxt = "n", 
#     bty = "n", xlab = "", las = 2,
#     main = datSet, lwd = 4)

#abline(h = seq(0,50,10), col = 8, lty = 3)

#par(mar = c(2,5,1,1))

#plot(0,0, type = "n",
#     xlim = c(0,length(samp)),
#     ylim = c(-1.5,ncol(filter_rm.sort)),
#     xaxs = "i", yaxs = "i", axes = FALSE,
#     xlab = "", ylab = "")

#for(i in 1:ncol(filter_rm.sort)){
#  select_filt <- (1:nrow(filter_rm.sort))[filter_rm.sort[,i]]
#  if(length(select_filt) > 0){
#    rect(xleft = select_filt - 1,
#         xright = select_filt,
#         ybottom = i - 1,
#         ytop = i, 
#         col = "black", 
#         density = -1, border = FALSE)
#  }
#}

#sum.cols <- RColorBrewer::brewer.pal(7,"GnBu")

#failedSamp <- rowSums(filter_rm.sort)
#failedSamp.no <- which(failedSamp > 0)
#failedSamp.col <- sum.cols[failedSamp[failedSamp>0]]

#rect(xleft = failedSamp.no - 1,
#     xright = failedSamp.no,
#     ytop = -.5,
#     ybottom = -1.5,
#     col = failedSamp.col, density = -1,
#     border = FALSE)

#mtext(side = 2, at = 1:ncol(filter_rm.sort) - .5,
#      text = colnames(filter_rm.sort), las = 2, line = .5)
#mtext(side = 2, at = -1,
#      text = "Total", las = 2, line = .5)
#abline( v= 0:1000, col= 8, lty = 3)


#plot(0,0, type = "n",
#     xlim = c(0,length(samp)),
#     ylim = c(0,3),
#     xaxs = "i", yaxs = "i", axes = FALSE,
#     xlab = "", ylab = "")
#

#x.vals <- 1:nrow(filter_rm.sort)
#names(x.vals) <- rownames(filter_rm.sort)
#trio_samps$x <- x.vals[trio_samps$Sample]

#trio_samps2 <- trio_samps[complete.cases(trio_samps),]

#rect(xleft = trio_samps2$x[trio_samps2$Group == "Chernobyl"] - 1,
#     xright = trio_samps2$x[trio_samps2$Group == "Chernobyl"],
#     ybottom = 2,ytop = 3, density = -1, col = 1)

#cols <- rep(1, sum(trio_samps2$Group == "Dog"))
#cols[trio_samps2$Source[trio_samps2$Group == "Dog"] == "Public"] <- 2
#cols[trio_samps2$Source[trio_samps2$Group == "Dog"] == "Dog_10k"] <- 4

#rect(xleft = trio_samps2$x[trio_samps2$Group == "Dog"] - 1,
#     xright = trio_samps2$x[trio_samps2$Group == "Dog"],
#     ybottom = 1,ytop = 2, density = -1, col = cols,
#     border = FALSE)

#cols <- rep(1, sum(trio_samps2$Group == "Wolf"))
#cols[trio_samps2$Source[trio_samps2$Group == "Wolf"] == "Public"] <- 2
#cols[trio_samps2$Source[trio_samps2$Group == "Wolf"] == "Dog_10k"] <- 4

#rect(xleft = trio_samps2$x[trio_samps2$Group == "Wolf"] - 1,
#     xright = trio_samps2$x[trio_samps2$Group == "Wolf"],
#     ybottom = 0,ytop = 1, density = -1, col = cols,
#     border = FALSE)

#mtext(side = 2, at = (3:1) - .5,
#      text = c("Chernobyl","Dog","Wolf"), las = 2, line = .5)

#legend("topright", fill = c(1,2,4), legend = c("In_house","Public","Dog_10k"), bty = "n",title = "Source", cex =.5)

#mtext(text = trio_samps2$Sample,side = 1,at = trio_samps2$x-.5,las = 2,cex = .4, line = .5)

#abline( v= 0:1000, col= 8, lty = 3)

#dev.off()





#write.table(data.frame(trio_samps,filter_rm.sort[trio_samps$Sample,], 
#                       PASS = rowSums(filter_rm.sort[trio_samps$Sample,]) == 0),
#            file = "~/Desktop/oneDrive/project_data/joint_coding_regions/trios/trio_samples_2216g_failed_filters.txt",
#            sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
#}



