
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("BiocGenerics")
BiocManager::install("GEOquery")

#library(BiocGenerics)
library(GEOquery)
#library(limma)
library(limma)

#check site with original data
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38739

#upload of data from NCBI website
gset <- getGEO("GSE38739", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6096", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

#make simple table with all expressions and names
my_table<-exprs(gset)

colnames(my_table)<-gset$title
genames<-gset@featureData@data$ILMN_Gene
row.names(my_table)<-genames

# write data to file
write.csv(my_table,file="GSE42414_fromGEO_R.csv")
write.table(my_table,file="GSE42414_fromGEO_R.tsv", sep="\t")

#some plots
title="set GSE42414"
boxplot(my_table, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2) 
hist_table<-hist(my_table, xlab="bins", breaks=50)
lines(hist_table$mids, hist_table$counts, "l", col="red", lwd=2)
plot(hist_table$mids, hist_table$counts, "l")

#for a fun make a subset of a data with higher expression only
mytable<-subset(my_table, rowMeans(my_table)>6 & rowMeans(my_table)<10)
dim(mytable)

#try to make correct colors, make and use shorter names
#get first part of the column names
library(dplyr)
library(purrr)
split<-map(strsplit(colnames(my_table), "_"),1) %>% unlist
split
table(split)
library(limma)
mds<-plotMDS(mytable, 
             main="Unknown expression data",
             labels=colnames(my_table),
             # pch=19,
             # col=split
             col=c(rep(1,5), rep(2,5), rep(3,5), 
                   rep(4,4),rep(5,5),rep(6,5),
                   rep(7,3),rep(8,3),rep(9,3), rep(10,3))
)

#much easier with colors in ggplot
library(ggplot2)
df<-data.frame(mds$x, mds$y, split)
ggplot(df, aes(x=mds.x, y=mds.y, col=split))+
  geom_point(cex=3)+
  ggtitle("Unknown expression data")

