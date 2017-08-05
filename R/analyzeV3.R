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

#________________________________ANALYZE_______________________________________
if (TRUE){
  #Variance & Fst
  varMatrix.alt = varMatrix.var = varMatrix.tot = 
    pBarMatrix.alt = pBarMatrix.var = pBarMatrix.tot =
    pFstMatrix.alt.aor = pFstMatrix.alt.roa =
    pFstMatrix.var = pFstMatrix.tot = 
    matrix(data=0, nrow=length(populations), ncol=length(populations))
  
  pVector = pVector.median = c()
  pFstMatrix.alt.median = matrix(nrow=0, ncol=length(populations))
}

for (r in 1:nrow(mer)){
  print(r)
  pVector = pVector.median = c()
  row.alt = as.numeric(altTable[r,])
  row.var = as.numeric(varTable[r,])
  row.tot = as.numeric(totTable[r,])
  
  for (c in 1:length(populations)){
    if (is.na(row.alt[c])) row.alt[c] = 0
    if (is.na(row.var[c])) row.var[c] = 0
    if (is.na(row.tot[c])) row.tot[c] = 0
  }
  P.alt = row.alt/(sizes*2)
  P.var = row.var/sizes
  P.tot = row.tot/sizes
  
  #w = sizes/sum(sizes)
  for (i in 1:length(populations)){
    for (j in 1:length(populations)){
      # Compute variance.
      size = sizes[i] + sizes[j]
      W.i = sizes[i] / size
      W.j = sizes[j] / size
      var.alt.ij = var(c(P.alt[i], P.alt[j])) / 2
      var.var.ij = var(c(P.var[i], P.var[j])) / 2
      var.tot.ij = var(c(P.tot[i], P.tot[j])) / 2
      
      varMatrix.alt[i,j] = varMatrix.alt[i,j] + var.alt.ij
      varMatrix.var[i,j] = varMatrix.var[i,j] + var.var.ij
      varMatrix.tot[i,j] = varMatrix.tot[i,j] + var.tot.ij
      
      # Compute Fst.
      pBar.alt.ij = sum(W.i * P.alt[i], W.j * P.alt[j])
      pBar.var.ij = sum(W.i * P.var[i], W.j * P.var[j])
      pBar.tot.ij = sum(W.i * P.tot[i], W.j * P.tot[j])

      pBarMatrix.alt[i,j] = pBarMatrix.alt[i,j] + pBar.alt.ij
      pBarMatrix.var[i,j] = pBarMatrix.var[i,j] + pBar.var.ij
      pBarMatrix.tot[i,j] = pBarMatrix.tot[i,j] + pBar.tot.ij
      
      Fst.ij = var.alt.ij / (pBar.alt.ij * (1 - pBar.alt.ij))
      #Fst.ij = ((P[i] - P[j]) ^ 2) / (2 * P.avg * (1 - P.avg))
      if (is.nan(Fst.ij)) Fst.ij = 0

      pFstMatrix.alt.aor[i,j] = pFstMatrix.alt.aor[i,j] + Fst.ij
      
      if (i != j) pVector = c(pVector, Fst.ij)
    }

    pVector[i] = mean(pVector[i:(i+24)])
    pVector.median[i] = median(pVector[i:(i+24)])
    pVector = pVector[1:i]
  }
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
  
  pBarMatrix.alt = pBarMatrix.alt / numRows
  pBarMatrix.var = pBarMatrix.var / numRows
  pBarMatrix.tot = pBarMatrix.tot / numRows
  
  pFstMatrix.alt.roa = varMatrix.alt / (pBarMatrix.alt * (1 - pBarMatrix.alt))
  pFstMatrix.alt.aor = pFstMatrix.alt.aor / numRows
  pFstMatrix.var = varMatrix.var / (pBarMatrix.var * (1 - pBarMatrix.var))
  pFstMatrix.tot = varMatrix.tot / (pBarMatrix.tot * (1 - pBarMatrix.tot))
  
  
  colnames(varMatrix.alt) = rownames(varMatrix.alt) = populations
  colnames(varMatrix.var) = rownames(varMatrix.var) = populations
  colnames(varMatrix.tot) = rownames(varMatrix.tot) = populations

  colnames(pFstMatrix.alt.roa) = rownames(pFstMatrix.alt.roa) = populations
  colnames(pFstMatrix.alt.aor) = rownames(pFstMatrix.alt.aor) = populations
  colnames(pFstMatrix.var) = rownames(pFstMatrix.var) = populations
  colnames(pFstMatrix.tot) = rownames(pFstMatrix.tot) = populations
  
  colnames(pFstMatrix.alt.avg) = colnames(pFstMatrix.alt.median) = populations
  rownames(pFstMatrix.alt.avg) = rownames(pFstMatrix.alt.median) = 
    paste(mer$CHR, ":", mer$POS, sep = "")
  
}

#____________________________________PLOT_______________________________________
if (TRUE){
  # Define Matrices and vectors.
  myMatricesLabels = c("Genetic distance based on variance of\nallele frequencies",
                       "Genetic distance based on variance of\nvariation frequencies",
                       "Genetic distance based on variance of\ntotal variation frequencies",
                       "Genetic distance based on Fst of\nallele frequencies",
                       "Genetic distance based on Fst of\nvariation frequencies",
                       "Genetic distance based on Fst of\ntotal variation frequencies",
                       "Fst of the frequency of a specific allele at a specific position in\neach population against all other populations at each position.")
  
}
# Plot.....
# Pairwise variance plots.
M2 = max(varMatrix.alt, varMatrix.var, varMatrix.tot)
myImagePlot(varMatrix.alt, min=0, max=M2)
title(main = myMatricesLabels[1], cex.main = 0.9)
myImagePlot(varMatrix.var, min=0, max=M2)
title(main = myMatricesLabels[2], cex.main = 0.9)
myImagePlot(varMatrix.tot, min=0, max=M2)
title(main = myMatricesLabels[3], cex.main = 0.9)

M = max(pFstMatrix.alt.aor, pFstMatrix.alt.roa, pFstMatrix.var, pFstMatrix.tot)
M = max(pFstMatrix.alt.aor, pFstMatrix.alt.roa)
# Pairwise Fst plots (average).
myImagePlot(pFstMatrix.alt.aor, min=0, max=M)
myImagePlot(pFstMatrix.alt.roa, min=0, max=M)
myImagePlot(pFstMatrix.alt.roa)
title(main = myMatricesLabels[4], cex.main = 0.9)
M = max(pFstMatrix.alt.roa, pFstMatrix.var, pFstMatrix.tot)
myImagePlot(pFstMatrix.var, min=0, max=M)
title(main = myMatricesLabels[5], cex.main = 0.9)
myImagePlot(pFstMatrix.tot, min=0, max=M)
title(main = myMatricesLabels[6], cex.main = 0.9)

# Fst plot (one subpopulation against the rest of the population) 
myImagePlot(pFstMatrix.alt.avg, border = FALSE)
myImagePlot(pFstMatrix.alt.median, border = FALSE)
filtered = pFstMatrix.alt.median[apply(pFstMatrix.alt.median[,1:26], MARGIN = 1, 
                                       function(x) any(x >= 0.35)),]
myImagePlot(filtered, border = TRUE, lwd=0.7)
chi.dist = numeric(26)
for (i in 1:length(populations)){
  for (j in 1: length(populations))
    if (i != j){ 
      chi.dist[i] = chi.dist[i] + ((r[i] - r[j]) ^ 2)
  }
}
chi.dist = chi.dist / 2
names(chi.dist) = populations

row=23
filtered[row, filtered[row,] >= 0.35]
M = max(pFstMatrix.alt.median)
chr = 1
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
clusDendro = dendrapply(hcp, colLab)
# make plot
plot(clusDendro)
