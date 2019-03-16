# qc_histogram for examining ATAC-seq insertion size

library(ggplot2)
library(plyr)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

# filename
hist_data <- read.delim(args[1])
# for re-binning
metric_value <- args[2]
# bin weights/reads
metric_count <- args[3]
# re-binning size
bin_size <- as.numeric(args[4])
# for calculating moving averages
smooth_size <- as.numeric(args[5])
# x axis title in plot
x_name <- args[6]
# plot output name (pdf)
out_name <- args[7]


# data must be sorted
hist_data = hist_data[order(hist_data[,metric_value]),]
# find re-bins
hist_data$qc_bin = round_any(hist_data[,metric_value], bin_size)
hist_data$qc_cum_quantile = cumsum(hist_data[,metric_count]) / sum(hist_data[,metric_count])

# only plot up to 99% of the data
outlier_threshold = head(subset(hist_data, qc_cum_quantile > 0.99),1)[,metric_value]
# find median and mean
median_size = head(subset(hist_data, qc_cum_quantile >= 0.5),1)[,metric_value]
mean_size = round(sum(hist_data[,metric_value] * (hist_data[,metric_count] / sum(hist_data[,metric_count]))))

# calculate re-bins
bin_data = ddply(hist_data, ~qc_bin, function(d) data.frame(qc_count = sum(d[,metric_count])))

moving_average = function(x, n = smooth_size){filter(x, rep(1/n, n), sides=2)}

# do not smooth at the number of reads when below 4 x the mean insert size or the max
mean_x_4 = round(4 * mean_size)
if(mean_x_4 > max(hist_data[,metric_value])) {
  threshold_metric_count = 1
} else {
  threshold_metric_value = mean_x_4 
  threshold_metric_count = head(hist_data[hist_data[,metric_value] >= threshold_metric_value, metric_count], 1) # should be small
}
thresholded_hist_data = hist_data[hist_data[,metric_count] > threshold_metric_count,]

# smooth only across bins with adequate coverage
if (smooth_size <= nrow(thresholded_hist_data)) {
  thresholded_hist_data$qc_smooth = moving_average(thresholded_hist_data[,metric_count])
} else {
  thresholded_hist_data$qc_smooth = thresholded_hist_data[,metric_count]
}
# normalize cdf to the extend of the y-axis
cdf_factor = max(thresholded_hist_data$qc_smooth,na.rm=TRUE)
hist_data$qc_norm_cum_quantile = hist_data$qc_cum_quantile * cdf_factor

#qc_color = "#3A5795"
qc_color <- "red"
window_size <- 600

pdf(out_name,width=4.2,height=3.3)

p1 <- (ggplot(bin_data, aes(x = qc_bin, y = qc_count)) 
  + theme_bw() 
  + theme(panel.grid.major=element_blank(), panel.grid.minor= element_blank()) 
  #+ geom_bar(stat="identity",fill="lightgray") 
  + geom_line(data = thresholded_hist_data, aes_string(x = metric_value, y = "qc_smooth"),color=qc_color,alpha=1) 
  #+ geom_line(data = hist_data, aes_string(x = metric_value, y = "qc_norm_cum_quantile"), color=qc_color,alpha=0.5) 
  #+ geom_vline(xintercept=mean_size, color=qc_color,alpha=0.5) 
  + xlab(x_name) 
  + ylab("Count") 
  + scale_y_continuous() 
  + scale_x_continuous(limits = c(0, window_size), 
                       breaks = seq(0, window_size, window_size/4), 
                       expand = c(0,0))
  #+ scale_y_log10()
  #+ scale_y_continuous(sec.axis = sec_axis(~./cdf_factor, name="CDF")) 
  + coord_cartesian(xlim=c(0, outlier_threshold))
)

p2 <- (ggplot(bin_data, aes(x = qc_bin, y = qc_count)) 
       + theme_bw() 
       + theme(panel.grid.major=element_blank(), panel.grid.minor= element_blank()) 
       #+ geom_bar(stat="identity",fill="lightgray") 
       + geom_line(data = thresholded_hist_data, aes_string(x = metric_value, y = "qc_smooth"),color=qc_color,alpha=1) 
       #+ geom_line(data = hist_data, aes_string(x = metric_value, y = "qc_norm_cum_quantile"), color=qc_color,alpha=0.5) 
       #+ geom_vline(xintercept=mean_size, color=qc_color,alpha=0.5) 
       + xlab(x_name) 
       + ylab("Count") 
       #+ scale_y_continuous() 
       + scale_y_log10()
       + scale_x_continuous(limits = c(0, window_size), 
                            breaks = seq(0, window_size, window_size/4), 
                            expand = c(0,0))
       #+ scale_y_continuous(sec.axis = sec_axis(~./cdf_factor, name="CDF")) 
       + coord_cartesian(xlim=c(0, outlier_threshold))
)

#grid.arrange(arrangeGrob(p1, p2, nrow=2), heights=c(4,1), widths=c(1,1))
p1
p2
dev.off()


