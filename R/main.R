# main.R
# A script to anylaze the VCF data and plot the results.
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
# Vectors describing the superpopulations.
superpopulations <- c("African", "American", "European", "East Asian", "South Asian")
# Start and end indices of each superpopulation in the subpopulations list.
super.start <- c(1, 8, 12, 17, 22)
super.end <- c(7, 11, 16, 21, 26)
# Make a list.
super = list("superpopulations" = superpopulations, 
             "super.s" = super.start, 
             "super.e" = super.end)

#______________________________EXTRACT_VARIATIONS_______________________________
getVariations(populations)

#______________________________MERGE_AND_ANALYZE________________________________
# Merge tables.
# 'tables': list of 4 tables: $altTable, $varTable, $totTable, and $coordinates.
tables = mergeData(populations, alt = TRUE, var = TRUE, tot = TRUE)
# Analyze.
# 'results': list of 12 tables.
results = analyzeData(populations, sizes, super, tables)

# Define data to be plotted...
# Coordinates needed for plotting pFstMatrix.alt.median by chromosome number.
coordinates = tables$coordinates

# Results matrices.
varMatrix.alt = results$varMatrix.alt
varMatrix.var = results$varMatrix.var
varMatrix.tot = results$varMatrix.tot

pFstMatrix.alt.roa = results$pFstMatrix.alt.roa
pFstMatrix.alt.aor = results$pFstMatrix.alt.aor
pFstMatrix.var = results$pFstMatrix.var
pFstMatrix.tot = results$pFstMatrix.tot

pFstMatrix.alt.median = results$pFstMatrix.alt.median

pFstMatrix.super = results$pFstMatrix.super

#____________________________________PLOT_______________________________________
# Define the labels for Matrices.
myLabels = c("Genetic distance based on variance of\nSNP frequencies",
             "Genetic distance based on variance of\nvariation frequencies",
             "Genetic distance based on variance of\ntotal-variation frequencies",
             "Genetic distance based on Fst of\nSNP frequencies",
             "Genetic distance based on Fst of\nvariation frequencies",
             "Genetic distance based on Fst of\ntotal-variation frequencies",
             "Fst estimates per SNP",
             "SNPs with high Fst values (>0.35)")

# Plot......
#______________________________Pairwise variance plots__________________________
varMax = max(varMatrix.alt, varMatrix.var, varMatrix.tot)

myImagePlot(varMatrix.alt, min=0, max=varMax, cex=0.8, title = myLabels[1])
myImagePlot(varMatrix.var, min=0, max=varMax, cex=0.8, title = myLabels[2])
myImagePlot(varMatrix.tot, min=0, max=varMax, cex=0.8, title = myLabels[3])

#______________________________Pairwise Fst plots_______________________________
# Plot pFstMatrix.alt.roa by itself.
myImagePlot(pFstMatrix.alt.roa, border = TRUE, sub = FALSE, 
            cex=0.8, title = myLabels[4])

# Plot pFstMatrix.alt.roa vs. pFstMatrix.alt.aor
fstMax_aoroa = max(pFstMatrix.alt.aor, pFstMatrix.alt.roa)
myImagePlot(pFstMatrix.alt.aor, min = 0, max = fstMax_aoroa,
            cex=0.8, title = myLabels[4])
myImagePlot(pFstMatrix.alt.roa, min = 0, max = fstMax_aoroa,
            cex=0.8, title = myLabels[4])

# Plot pFstMatrix.alt.roa vs. pFstMatrix.var vs. pFstMatrix.tot
fstMax_altvartot = max(pFstMatrix.alt.roa, pFstMatrix.var, pFstMatrix.tot)
myImagePlot(pFstMatrix.alt.roa, min = 0, max = fstMax_altvartot, 
            cex=0.8, title = myLabels[1])
myImagePlot(pFstMatrix.var, min = 0, max = fstMax_altvartot, 
            cex=0.8, title = myLabels[5])
myImagePlot(pFstMatrix.tot, min = 0, max = fstMax_altvartot, 
            cex=0.8, title = myLabels[6])

#______________________________pFstMatrix.super_________________________________
myImagePlot(pFstMatrix.super, cex=1, title = myLabels[4])

#_____________________________pFstMatrix.alt.median_____________________________
myImagePlot(pFstMatrix.alt.median, border = FALSE, cex=0.8, title = myLabels[7])

# Plot pFstMatrix.alt.median by chromosome number.
# (change the value of chr to be the chromosome number.)
M = max(pFstMatrix.alt.median)
chr = 22
myImagePlot(pFstMatrix.alt.median[coordinates$CHR %in% chr, ], 
            lwd=0.2, min = 0, max = M, cex = 0.7, title = myLabels[7])

# Filter rows containing at least one Fst value > 0.35
filtered = pFstMatrix.alt.median[apply(pFstMatrix.alt.median[,1:26], MARGIN = 1, 
                                       function(x) any(x >= 0.35)),]
myImagePlot(filtered, border = TRUE, cex = 0.8, title = myLabels[8])


#________________________________Dendrogram_____________________________________
par(mar=c(5,5,1,1))

subDendrogram = getDendrogram(pFstMatrix.alt.roa, 6)
plot(subDendrogram, ylab = "Fst")

superDendrogram = getDendrogram(pFstMatrix.super, 5)
plot(superDendrogram, ylab = "Fst")
