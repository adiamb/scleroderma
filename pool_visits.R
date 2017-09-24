args = commandArgs(trailingOnly=TRUE)
files = args

package_path = '/path/to/packages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)

data = data.frame()
for (file in files) {
  week = strsplit(file,'_')[[1]][2]
  data1$WEEK = week
  data1$WEEK_CORRECTED = convert_nominal_time_to_precise_time(week) # convert_nominal_time_to_precise_time: some function
  data = rbind(data,data1)
}

data = data[order(data$WEEK),]
outfile = gsub(paste("_",week,sep=""),"",file)
outfile = gsub(".tab","_visits-pooled.tab",outfile)

writeChangeoDb(data,outfile)