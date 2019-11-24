# Rscript calculate_counts_from_depth_and_vaf_for_one_sample.r infile outfile

# install.packages("tvd")

options(scipen = 999)
options(width=150)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))    # For plotting only
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tvd))        # For plotting only

args = commandArgs(trailingOnly=TRUE) # for production
#args=c("./AABSG_LP1001138-NTP_H10_ST-E00110_HYH3JCCXX_3.loh.nb2.rds", "AABSG_LP1001138-NTP_H10_ST-E00110_HYH3JCCXX_3_nb2_counts.tsv") # for testing

infile=args[1]
outfile=args[2]

return_data_from_plotWindowFit = function(data, fit)
{
    fit$kf1 = fit$fit.f*fit$fit.k1/2 + (1-fit$fit.f)/2
    fit$kf2 = fit$fit.f*fit$fit.k2/2 + (1-fit$fit.f)/2
    fit$d = log2(fit$kf1 + fit$kf2)
    fit$v1 = fit$kf1 / (fit$kf1 + fit$kf2)
    fit$v2 = fit$kf2 / (fit$kf1 + fit$kf2)
    fit$colus = fit$fit.llik / (fit$end_index - fit$start_index + 1)
    fit$col = pmax(0, pmin(1, (fit$colus - (-10)) / 5))
    fit$coldisc = as.integer(ceiling(fit$col * 100))
    fit$coldisc[fit$coldisc == 0] = 1
    pal = rainbow(100, end = 0.33)

    chrom_boundaries = which(data$chrom[-1] != data$chrom[-nrow(data)])
    chrom_midpoints = (c(0, chrom_boundaries) + c(chrom_boundaries, nrow(data))) / 2
    names(chrom_midpoints) = data$chrom[!duplicated(data$chrom)]

    #par(mfrow = c(2, 1))
    #plot(log2(data$dp) - log2(data$pois.lambda), pch = ".", col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "Depth anomaly (log2)", main = "Depth", xaxt = "n") # <== black cloud of dots
    #axis(1, at = chrom_midpoints, labels = names(chrom_midpoints)) # <== x-axis chromosome numbers
    #lines(tvd1d(log2(data$dp) - log2(data$pois.lambda), lambda = 5), col = "red", lwd = 2) # <== the red horizontal bits of lines under the green horizontal line
    #abline(v = chrom_boundaries, col = "blue", lwd = 2) # <== dark blue vertical lines
    #abline(v = fit$start_index, col = "red", lty = "dotted") # <== red vertical dotted lines
    #segments(fit$start_index, fit$d, fit$end_index, fit$d, col = pal[fit$coldisc], lwd = 3) # <== this is the middle green line

    #plot(data$ad / jitter(data$dp, amount = 0.5), pch = ".", ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "VAF", main = "Allele frequency", xaxt = "n") # <== black cloud of dots
    #axis(1, at = chrom_midpoints, labels = names(chrom_midpoints)) # <== x-axis chromosome numbers
    #abline(v = chrom_boundaries, col = "blue", lwd = 2) # <== dark blue vertical lines
    #abline(v = fit$start_index, col = "red", lty = "dotted") # <== red vertical dotted lines
    #segments(fit$start_index, fit$v1, fit$end_index, fit$v1, col = pal[fit$coldisc], lwd = 3) # <== this is the middle green horizontal line and top little green line bits
    #segments(fit$start_index, fit$v2, fit$end_index, fit$v2, col = pal[fit$coldisc], lwd = 3) # <== this is the middle green horizontal line and bottom little green line bits
    #par(mfrow = c(1, 1))

    start_chrom = data[ fit$start_index, c("chrom") ]
    end_chrom = data[ fit$end_index, c("chrom") ]
    green_lines = cbind( start_chrom, end_chrom, fit$start_index, fit$end_index, fit$d, fit$v1, fit$v2, fit$fit.f )
    colnames(green_lines) = c( "start_chrom", "end_chrom", "start_index", "end_index", "depth", "vaf1", "vaf2", "measure_of_subclonality_from_model_fitting" )
    #red_lines = tvd1d(log2(data$dp) - log2(data$pois.lambda), lambda = 5) # has the same number of elements as data has rows = 713365
    return_data = green_lines
}

nb2=readRDS( infile, refhook = NULL)
data = nb2$data
fit = nb2$fit
green_lines = return_data_from_plotWindowFit(data, fit$fit)
write.table(green_lines, file=outfile, row.names=FALSE, na="", col.names=TRUE, sep="\t", quote=FALSE)


