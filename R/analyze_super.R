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
  varMatrix.alt = pBarMatrix.alt = pFstMatrix.alt = pFstMatrix.alt.aor =
    matrix(data=0, nrow=length(super), ncol=length(super))
}

for (r in 1:nrow(mer)){
  cat(r, "...")
  row.alt = as.numeric(altTable.super[r,])

  for (c in 1:length(super)){
    if (is.na(row.alt[c])) row.alt[c] = 0
  }
  P.alt = row.alt / (sizes.super * 2)
  for (i in 1:length(super)){
    for (j in 1:length(super)){
      size = sizes.super[i] + sizes.super[j]
      W.i = sizes.super[i] / size
      W.j = sizes.super[j] / size
      # Compute variance.
      var.alt.ij = var(c(P.alt[i], P.alt[j])) / 2
      varMatrix.alt[i,j] = varMatrix.alt[i,j] + var.alt.ij
      # Compute pBar
      pBar.alt.ij = sum(W.i * P.alt[i], W.j * P.alt[j])
      pBarMatrix.alt[i,j] = pBarMatrix.alt[i,j] + pBar.alt.ij
      
      Fst.ij = var.alt.ij / (pBar.alt.ij * (1 - pBar.alt.ij))
      
      if (is.nan(Fst.ij)) Fst.ij = 0
      pFstMatrix.alt.aor[i,j] = pFstMatrix.alt.aor[i,j] + Fst.ij
    }
  }
}
#__________________________________________________
#varListbyCHR = varListbyCHR[2:length(varListbyCHR)]
if (TRUE){
  numRows = nrow(mer)
  varMatrix.alt = varMatrix.alt / numRows
  pBarMatrix.alt = pBarMatrix.alt / numRows
  pFstMatrix.alt.aor = pFstMatrix.alt.aor / numRows
  
  pFstMatrix.alt = varMatrix.alt / (pBarMatrix.alt * (1 - pBarMatrix.alt))
  colnames(pFstMatrix.alt) = rownames(pFstMatrix.alt) = super
}
#____________________________________PLOT_______________________________________
# Plot.....
# Pairwise Fst plots (average).
myImagePlot(pFstMatrix.alt)
title(main = myMatricesLabels[4], cex.main = 0.9)
myImagePlot(pFstMatrix.alt.aor)
