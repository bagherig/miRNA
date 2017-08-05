#Merge
path = "/home/pinterbagz/Desktop/Research/VCF/"
fileName = "/results.tsv"
populations <- c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI",  # [1:7] AFRICAN
                 "CLM", "MXL", "PEL", "PUR",                       # [8:11] AMERICAN
                 "CEU", "FIN", "GBR", "IBS", "TSI",                # [12:16] EUROPEAN
                 "CDX", "CHB", "CHS", "JPT", "KHV",                # [17,21] EAST ASIAN
                 "BEB", "GIH", "ITU", "PJL", "STU")                # [22:26]SOUTH ASIAN
sizes <- c(95, 61, 99, 113, 99, 85, 108, 
           94, 64, 85, 104, 
           99, 99, 91, 107, 107, 
           93, 103, 105, 104, 99, 
           86, 103, 102, 96, 102)
super <- c("African", "American", "European", "East Asian", "South Asian")
super.s <- c(1, 8, 12, 17, 22)
super.e <- c(7, 11, 16, 21, 26)
sizes.super = c(sum(sizes[1:7]), sum(sizes[8:11]), sum(sizes[12:16]), 
                sum(sizes[17:21]), sum(sizes[22:26]))

numChr = 22

#________________________________Merge Data_______________________________________
mer <- NULL
for (pop in populations){
  cat ("Merging", pop, "... ")
  table <- read.table(paste(path, pop, fileName, sep = ""), sep = "\t", 
                      header = TRUE, stringsAsFactors = FALSE)
  table <- table[table$totFREQ >= 5, c("CHR", "POS", "ALT", "altFREQ", "varFREQ", "totFREQ")]
  if (is.null(mer)){
    mer <- table
  } else {
    mer <- merge(mer, table, all = TRUE, 
                 by = c("CHR", "POS", "ALT"))
  }
  mer = mer[!(duplicated(mer)),]
}

altTable = mer[, seq(4, ncol(mer), 3)]
varTable = mer[, seq(5, ncol(mer), 3)]
totTable = mer[, c(1, 2, seq(6, ncol(mer), 3))]
totTable = totTable[!(duplicated(totTable)), 3:ncol(totTable)]
colnames(altTable) = colnames(varTable) = colnames(totTable) = populations

altTable.super = matrix(ncol = 5, nrow = nrow(altTable))
colnames(altTable.super) = super
for (i in 1:26){
  for (j in 1:nrow(altTable)){
    if (is.na(altTable[j,i])){ altTable[j,i]=0}
  }
}
for (i in 1:5){
  altTable.super[,i] = apply(altTable[,super.s[i]:super.e[i]], 1, sum)
}
#________________________________ANALYZE_______________________________________
if (TRUE){
  varMatrix.alt = varMatrix.var = varMatrix.tot = 
    pFstMatrix.alt = pFstMatrix.var = pFstMatrix.tot = 
    matrix(data=0, nrow=length(super), ncol=length(super))
  
  sVector = pVector = pVector.median = c()
  pFstMatrix.alt.avg = pFstMatrix.alt.median = matrix(nrow=0, ncol=length(super))
}

for (r in 1:nrow(mer)){
  cat(r, "...")
  sVector = pVector = c()
  row.alt = as.numeric(altTable.super[r,])
  row.var = as.numeric(varTable[r,])
  row.tot = as.numeric(totTable[r,])
  for (c in 1:length(super)){
    if (is.na(row.alt[c])) row.alt[c] = 0
    if (is.na(row.var[c])) row.var[c] = 0
    if (is.na(row.tot[c])) row.tot[c] = 0
  }
  freq.alt = row.alt/sizes.super/2
  freq.var = row.var/sizes.super
  freq.tot = row.tot/sizes.super
  
  for (i in 1:length(super)){
    for (j in 1:length(super)){
      # Compute variance.
      var.alt.ij = var(c(freq.alt[i], freq.alt[j]))
      var.var.ij = var(c(freq.var[i], freq.var[j]))
      var.tot.ij = var(c(freq.tot[i], freq.tot[j]))
      
      varMatrix.alt[i,j] = varMatrix.alt[i,j] + var.alt.ij
      varMatrix.var[i,j] = varMatrix.var[i,j] + var.var.ij
      varMatrix.tot[i,j] = varMatrix.tot[i,j] + var.tot.ij
      
      # Compute Fst.
      #varSize.i = l[i] * sizes[i] / 100
      #varSize.j = l[j] * sizes[j] / 100
      #pBar.alt.ij = (row.alt[i] + row.alt[j]) / (sizes[i]*2 + sizes[j]*2)
      pBar.alt.ij = (freq.alt[i] + freq.alt[j]) / 2
      pBar.var.ij = (freq.var[i] + freq.var[j]) / 2
      pBar.tot.ij = (freq.tot[i] + freq.tot[j]) / 2
      
      #Fst.alt.ij = var.alt.ij / (pBar.alt.ij * (1 - pBar.alt.ij))
      Fst.alt.ij = var.alt.ij / pBar.alt.ij
      Fst.var.ij = var.var.ij / pBar.var.ij
      Fst.tot.ij = var.tot.ij / pBar.tot.ij
      if (is.nan(Fst.alt.ij)) Fst.alt.ij = 0
      if (is.nan(Fst.var.ij)) Fst.var.ij = 0
      if (is.nan(Fst.tot.ij)) Fst.tot.ij = 0
      
      pFstMatrix.alt[i,j] = pFstMatrix.alt[i,j] + Fst.alt.ij
      pFstMatrix.var[i,j] = pFstMatrix.var[i,j] + Fst.var.ij
      pFstMatrix.tot[i,j] = pFstMatrix.tot[i,j] + Fst.tot.ij
      
      if (i != j) pVector = c(pVector, Fst.alt.ij)
    }
    pVector[i] = mean(pVector[i:(i+3)])
    pVector.median[i] = median(pVector[i:(i+3)])
    pVector = pVector[1:i]
  }
  sFstMatrix.alt = rbind(sFstMatrix.alt, sVector)
  pFstMatrix.alt.avg = rbind(pFstMatrix.alt.avg, pVector)
  pFstMatrix.alt.median = rbind(pFstMatrix.alt.median, pVector.median)
  
}
#__________________________________________________
#varListbyCHR = varListbyCHR[2:length(varListbyCHR)]
if (TRUE){
  numRows = nrow(mer)
  varMatrix.alt = varMatrix.alt / numRows
  varMatrix.var = varMatrix.var / numRows
  varMatrix.tot = varMatrix.tot / numRows
  #chiMatrix = chiMatrix / nrow(mer)
  pFstMatrix.alt = pFstMatrix.alt / numRows
  pFstMatrix.var = pFstMatrix.var / numRows
  pFstMatrix.tot = pFstMatrix.tot / numRows
  
  FstBarPlot.alt = FstBarPlot.alt / numRows
  FstBarPlot.var = FstBarPlot.var / numRows
  FstBarPlot.tot = FstBarPlot.tot / numRows
  
  colnames(varMatrix.alt) = rownames(varMatrix.alt) = super
  colnames(varMatrix.var) = rownames(varMatrix.var) = super
  colnames(varMatrix.tot) = rownames(varMatrix.tot) = super
  #colnames(chiMatrix) = rownames(chiMatrix) = super
  colnames(pFstMatrix.alt) = rownames(pFstMatrix.alt) = super
  colnames(pFstMatrix.var) = rownames(pFstMatrix.var) = super
  colnames(pFstMatrix.tot) = rownames(pFstMatrix.tot) = super
  
  colnames(pFstMatrix.alt.avg) = colnames(pFstMatrix.alt.median) = super
  rownames(pFstMatrix.alt.avg) = rownames(pFstMatrix.alt.median) = 
    paste(mer$CHR, ":", mer$POS, sep = "")
}
#____________________________________PLOT_______________________________________
if (TRUE){
  # Define Colors.
  matColRamp = rgb(seq(1, 1, length=256),  # Red
                   seq(1, 0, length=256),  # Green
                   seq(1, 0, length=256))  # Blue
  vecColRamp = rgb(seq(0, 1, length=length(super)),  # Red
                   seq(0.5, 0, length=length(super)),  # Green
                   seq(1, 0, length=length(super)))  # Blue
  chrColRamp = (rgb(seq(0, 1, length=2),  # Red
                    seq(0.5, 0, length=2),  # Green
                    seq(1, 0, length=2)))  # Blue
  chrColors = c()
  for (chr in 1:22){
    if (chr%%2 == 1){
      chrColors = c(chrColors, rep(chrColRamp[1], sum(mer$CHR == chr)))
    } else {
      chrColors = c(chrColors, rep(chrColRamp[2], sum(mer$CHR == chr)))
    }
  }
  # Define Matrices and vectors.
  myVectors = c("var.alt", "fst.alt",
                "FstBarPlot.alt", "FstBarPlot.var", "FstBarPlot.tot")
  myMatrices = c("varMatrix.alt", "varMatrix.var", "varMatrix.tot", 
                 "pFstMatrix.alt", "pFstMatrix.var", "pFstMatrix.tot",
                 "sFstMatrix.alt", "pFstMatrix.alt.avg", "pFstMatrix.alt.median")
  myMatricesLabels = c("Genetic distance based on variance of\nallele frequencies",
                       "Genetic distance based on variance of\nvariation frequencies",
                       "Genetic distance based on variance of\ntotal variation frequencies",
                       "Genetic distance based on Fst of\nallele frequencies",
                       "Genetic distance based on Fst of\nvariation frequencies",
                       "Genetic distance based on Fst of\ntotal variation frequencies",
                       "Fst of the frequency of a specific allele at a specific position in\neach population against all other super at each position.")
  
}
# Plot.....
# Pairwise variance plots.
M = max(pFstMatrix.alt, pFstMatrix.var, pFstMatrix.tot)
M2 = max(varMatrix.alt, varMatrix.var, varMatrix.tot)
myImagePlot(varMatrix.alt, min=0, max=M2)
title(main = myMatricesLabels[1], cex.main = 0.9)
myImagePlot(varMatrix.var, min=0, max=M2)
title(main = myMatricesLabels[2], cex.main = 0.9)
myImagePlot(varMatrix.tot, min=0, max=M2)
title(main = myMatricesLabels[3], cex.main = 0.9)

# Pairwise Fst plots (average).
myImagePlot(pFstMatrix.alt, min=0, max=M)
title(main = myMatricesLabels[4], cex.main = 0.9)
myImagePlot(pFstMatrix.var, min=0, max=M)
title(main = myMatricesLabels[5], cex.main = 0.9)
myImagePlot(pFstMatrix.tot, min=0, max=M)
title(main = myMatricesLabels[6], cex.main = 0.9)

# Fst plot (one subpopulation against the rest of the population) 
myImagePlot(pFstMatrix.alt.median, border = FALSE)
filtered = pFstMatrix.alt.median[apply(pFstMatrix.alt.median[,1:5], MARGIN = 1, 
                                       function(x) any(x >= 0.35)),]
myImagePlot(filtered, border = TRUE, lwd=0.7)
row=23
filtered[row, filtered[row,] >= 0.35]
M = max(pFstMatrix.alt.median)
chr = 1
myImagePlot(pFstMatrix.alt.median[mer$CHR %in% chr, ], 
            lwd=0.2, min = 0, max = M, cex = 0.5)

altTable[mer$CHR==4 & mer$POS==184851037, ]/sizes.super/2 ############################

title(main = myMatricesLabels[7], cex.main = 0.9)
