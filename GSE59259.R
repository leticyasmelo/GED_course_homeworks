if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("BiocGenerics")
BiocManager::install("GEOquery")

library(GEOquery)
library(limma)
library(limma)
library(GEOquery)
install.packages("pheatmap")
library(pheatmap)
install.packages("calibrate")
library(calibrate)

# Shows volcano plot and heatmap for 5 most DEGs.
VolcHeatPlot <- function (myfit, results, suffix = "adj") {
  
  # print (results["adj.P.Val"][1:10,])
  print (results[1:10,])
  write.csv(results["adj.P.Val"], sprintf("de-%s.csv", suffix))
  
  # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
  my_sample_col <- data.frame(sample = rep(c("normal", "tumour"), c(8,8)))
  row.names(my_sample_col) <- colnames(exprs(gset[[1]]))
  
  pheatmap(exprs(gset[[1]][rownames(results[1:5,])]),
           # annotation_row = ,
           annotation_col = my_sample_col,
           cutree_rows = 2,
           cutree_cols = 2,
           filename = sprintf("heatmap-%s.pdf", suffix)
  )
  
  # print (summary(results))
  # venResults <- decideTests(myfit)
  # vennDiagram(results)
  # vennCounts(venResults)
  
  # set parameters and draw the plot
  png(sprintf("VolcanoPlot%s.png", suffix))
  # mytitle <- paste ("GSE59259", '/', "selected samples ", sprintf("[%s]", suffix), sep ='')
  # X11(title = mytitle)
  # dev.new()
  
  # plotMA(myfit)
  volcanoplot(myfit)
  
  # taken from https://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
  
  # Add colored points: red if padj<0.01, orange of log2FC>2, green if both)
  with(results, points(logFC, -log10(P.Value), pch = 20, col = "green"))
  with(results[1:5,], points(
    logFC,-log10(P.Value),
    pch = 20,
    col = "red",
    cex = 2
  ))
  
  # from calibrate library
  with(results, textxy(
    logFC,-log10(P.Value),
    labs = rownames(results),
    cex = .8
  ))
  dev.off()
}


if (!exists("gset")) {
  gset <- getGEO("GSE59259", GSEMatrix = TRUE)
}

# making design matrix to fit linear model
type <- c(rep("0", 8), rep("1", 8))
type <- factor(type)
des <- model.matrix(~ 0 + type)
colnames(des) <- c("N", "H")

Data <- log(exprs(gset[[1]]) + 0.001)
fit <- lmFit(Data, design = des) ## !=

# boxplot and histogram about the data
title="data GSE59259"
boxplot(Data, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2) 
hist_table<-hist(Data, xlab="bins", breaks=100)
lines(hist_table$mids, hist_table$counts, "l", col="red", lwd=2)
plot(hist_table$mids, hist_table$counts, "l")

# making contrast matrix for differentially expressed genes
contrast.matrix <- makeContrasts(Diff = N - H, levels = des)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

t <- topTable(fit2, n = Inf)

# results <- t[t$adj.P.Val < 0.01 & abs(t$logFC) > 2,]

VolcHeatPlot (results = t[t$adj.P.Val < 0.01, ], myfit = fit2)
VolcHeatPlot (fit2, t[t$adj.P.Val < 0.01 & abs(t$logFC) > 2, ], "adjlogFC")

# analyzing the data we can observe differential expressed genes
