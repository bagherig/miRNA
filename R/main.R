# main.R
#
# A vector of subpopulation codes.
populations = c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", # [1:7] AFRICAN
                "CLM", "MXL", "PEL", "PUR",                    # [8:11] AMERICAN
                "CEU", "FIN", "GBR", "IBS", "TSI",            # [12:16] EUROPEAN
                "CDX", "CHB", "CHS", "JPT", "KHV",          # [17,21] EAST ASIAN
                "BEB", "GIH", "ITU", "PJL", "STU")          # [22:26]SOUTH ASIAN
# A vector of the sizes of subpopulations. The order must be the same as
# the "populations" vector.
sizes <- c(95, 61, 99, 113, 99, 85, 108, 
           94, 64, 85, 104, 
           99, 99, 91, 107, 107, 
           93, 103, 105, 104, 99, 
           86, 103, 102, 96, 102)

#______________________________MERGE_AND_ANALYZE________________________________

tables = mergeData(populations, alt = TRUE, var = TRUE, tot = TRUE)
results = analyzeData(populations, sizes, tables)

varMatrix.alt = results$varMatrix.alt
varMatrix.var = results$varMatrix.var
varMatrix.tot = results$varMatrix.tot

pFstMatrix.alt.roa = results$pFstMatrix.alt.roa
pFstMatrix.alt.aor = results$pFstMatrix.alt.aor
pFstMatrix.var = results$pFstMatrix.var
pFstMatrix.tot = results$pFstMatrix.tot

pFstMatrix.alt.median = results$pFstMatrix.alt.median


#____________________________________PLOT_______________________________________
# Define the labels for Matrices.
myMatricesLabels = c("Genetic distance based on variance of\nSNP frequencies",
                     "Genetic distance based on variance of\nvariation frequencies",
                     "Genetic distance based on variance of\ntotal-variation frequencies",
                     "Genetic distance based on Fst of\nSNP frequencies",
                     "Genetic distance based on Fst of\nvariation frequencies",
                     "Genetic distance based on Fst of\ntotal-variation frequencies",
                     "Fst estimates per SNP")

# Plot......
# Pairwise variance plots...
varMax = max(varMatrix.alt, varMatrix.var, varMatrix.tot)
myImagePlot(varMatrix.alt, min=0, max=varMax, cex=0.7)
title(main = myMatricesLabels[1], cex.main = 0.9)
myImagePlot(varMatrix.var, min=0, max=varMax)
title(main = myMatricesLabels[2], cex.main = 0.9)
myImagePlot(varMatrix.tot, min=0, max=varMax)
title(main = myMatricesLabels[3], cex.main = 0.9)

# Pairwise Fst plots...
# Plot pFstMatrix.alt.roa by itself.
myImagePlot(pFstMatrix.alt.roa)
fstMax_altvartot = max(pFstMatrix.alt.roa, pFstMatrix.var, pFstMatrix.tot)

# Plot pFstMatrix.alt.roa vs. pFstMatrix.alt.aor
fstMax_aoroa = max(pFstMatrix.alt.aor, pFstMatrix.alt.roa)
myImagePlot(pFstMatrix.alt.aor, min = 0, max = fstMax_aoroa)
myImagePlot(pFstMatrix.alt.roa, min = 0, max = fstMax_aoroa)

# Plot pFstMatrix.alt.roa vs. pFstMatrix.var vs. pFstMatrix.tot
myImagePlot(pFstMatrix.alt.roa)
title(main = myMatricesLabels[4], cex.main = 0.9)
myImagePlot(pFstMatrix.var, min=0, max=M)
title(main = myMatricesLabels[5], cex.main = 0.9)
myImagePlot(pFstMatrix.tot, min=0, max=M)
title(main = myMatricesLabels[6], cex.main = 0.9)

# Plot pFstMatrix.alt.median
myImagePlot(pFstMatrix.alt.avg, border = FALSE)
myImagePlot(pFstMatrix.alt.median, border = FALSE)
title(main = myMatricesLabels[7], cex.main = 0.9)
filtered = pFstMatrix.alt.median[apply(pFstMatrix.alt.median[,1:26], MARGIN = 1, 
                                       function(x) any(x >= 0.35)),]
myImagePlot(filtered, border = TRUE)

M = max(pFstMatrix.alt.median)
chr = 22
myImagePlot(pFstMatrix.alt.median[mer$CHR %in% chr, ], 
            lwd=0.2, min = 0, max = M, cex = 0.5)

altTable[mer$CHR==4 & mer$POS==184851037, ]/sizes/2 ############################

title(main = myMatricesLabels[7], cex.main = 0.9)

for (i in 1:26){
  for (j in i:26){
    p[i,j]=0
  }
}
myImagePlot(pFstMatrix.alt.roa, border = FALSE)
myImagePlot(pFstMatrix.alt.roa, lwd = 0.3)
myImagePlot(pFstMatrix.alt.roa, lwd = 0.3, sub=TRUE)
myImagePlot(p, border = FALSE)
myImagePlot(p, half = TRUE, lwd = 0.3)
myImagePlot(p, half = TRUE, lwd = 0.3, sub = TRUE)
heatmap(p)
f<-function(x){
  return(as.dist(x))
}
heatmap(pFstMatrix.alt.roa, col=ColorRamp, symm = TRUE, distfun = f)
heatmap(pFstMatrix.alt.roa, col=ColorRamp, distfun = f)

#Draw dendrogram
h=hclust(as.dist(pFstMatrix.alt.roa))
labelColors = c("blue", "red", "purple", "darkgreen", "brown")
clusMember = cutree(h, 6)

colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
# using dendrapply
hcd=as.dendrogram(h)
clusDendro = dendrapply(hcd, colLab)
# make plot
par(mar=c(4,5,1,0))
plot(clusDendro, ylab = "Fst")
