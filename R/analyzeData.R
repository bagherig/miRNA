#________________________________ANALYZE_______________________________________
#' Calculate the genetic distance between subpopulations based on variance and
#' fixation index. 
#'
#' @param populations A vector of subpopulation codes.
#' @param sizes A vector of sizes of each subpopulation. The order must be the
#'              same as "populations".
#' @param tables A list of tables with frequency values for each subpopulation.
#'
#' @return A list of matrices containing variance, pBar, and Fst values.
#'
#' @examples analyzeData(list(altTable, varTable, totTable))
#' 
#' @details This function was not designed for efficieny and speed. 
#' 
analyzeData <- function(populations = NULL, sizes = NULL, super = NULL, tables){
  # KEYWORDS: .alt -> SNP frequencies
  #           .var -> variation frequencies
  #           .tot -> total-variation frequencies
  #           .aor -> average of ratios
  #           .roa -> ratio of averages
  #           .super -> matrix for superpopulations
  #
  # Error checking:
  if (is.null(populations)){ 
    stop("a list of subpopulation codes must be specified!") }
  if (is.null(sizes)){ 
    stop("a list of subpopulation sizes must be specified!") }
  
  #___________________________ANALYZE SUBPOPULATIONS____________________________  
  # Initiate matrices and vectors to store the results.
  # Matrices for storing variances.
  varMatrix.alt = varMatrix.var = varMatrix.tot = 
    # Matrices for storing pBar.
    pBarMatrix.alt = pBarMatrix.var = pBarMatrix.tot = 
    # Matrices for storing Fst estimates. 
    # aor is "average of ratios", roa is "ratio of averages"
    pFstMatrix.alt.aor = pFstMatrix.alt.roa = 
    pFstMatrix.var = pFstMatrix.tot = 
    matrix(data=0, nrow=length(populations), ncol=length(populations))
  
  # Vector for storing the Fst estimates of each row for finding 
  # characteristic SNPs.
  pVector = c()
  # Matrix to store all pVectors.
  pFstMatrix.alt.median = matrix(nrow=0, ncol=length(populations))
  
  #_________________________________ANALYZE_____________________________________
  # Define input tables.
  coordinates = tables$coordinates
  altTable = varTable = totTable = NULL
  if ("altTable" %in% names(tables)){ 
    altTable = tables$altTable }
  if ("varTable" %in% names(tables)){ 
    varTable = tables$varTable }
  if ("totTable" %in% names(tables)){ 
    totTable = tables$totTable }
  
  numRows = max(nrow(altTable), nrow(varTable), nrow(totTable))
  # Go through every row of the table.
  for (r in 1:numRows){
    # Print the row number being analyzed.
    cat("Analyzing row", r, "\n")
    # Clear pVector for each row.
    pVector = c()
    
    if (!(is.null(altTable))){ 
      row.alt = as.numeric(altTable[r,]) # Read the current row.
      row.alt[is.na(row.alt)] = 0   # switch NA values to 0.
      P.alt = row.alt/(sizes * 2)     # calculate the frequencies of this SNP.
    }
    if (!(is.null(varTable))){ 
      row.var = as.numeric(varTable[r,]) # Read the current row.
      row.var[is.na(row.var)] = 0   # switch NA values to 0.
      P.var = row.var/sizes       # calculate the frequencies of this variation.
    }
    if (!(is.null(totTable))){ 
      row.tot = as.numeric(totTable[r,]) # Read the current row.
      row.tot[is.na(row.tot)] = 0   # switch NA values to 0.
      P.tot = row.tot/sizes  # calculate the frequencies of the total-variation.
    }
    
    # Doing pairwise comparison
    for (i in 1:length(populations)){
      for (j in 1:length(populations)){
        # Calculate the sum of the sizes of the two subpopulations.
        size = sizes[i] + sizes[j]
        # Calculate the proportion of size of each subpopulation.
        W.i = sizes[i] / size  
        W.j = sizes[j] / size
        
        if (!(is.null(altTable))){ 
          # Compute variance. Divide variance by 2 since the R function var()
          # divides by N-1 (N = 2 here) when calculating variance, but we want
          # to divide by N.
          var.alt.ij = var(c(P.alt[i], P.alt[j])) / 2
          # Add the calculated variance to the matrix of variances. We will
          # divide by the number of SNPs at the end.
          varMatrix.alt[i,j] = varMatrix.alt[i,j] + var.alt.ij
          
          # Compute pBar.
          pBar.alt.ij = sum(W.i * P.alt[i], W.j * P.alt[j])
          pBarMatrix.alt[i,j] = pBarMatrix.alt[i,j] + pBar.alt.ij
          
          # Computer Fst.
          Fst.ij = var.alt.ij / (pBar.alt.ij * (1 - pBar.alt.ij))
          # It is possible that both the variance and pBar are 0. In this case,
          # R returns NaN for 0/0, however, we want an Fst value of 0.
          if (is.nan(Fst.ij)) { Fst.ij = 0 }
          # Add to the matrix of Fst values. This matrix is for calculating 
          # genome-wide Fst by taking the average of ratios. We will divide by
          # the number of subpopulations at the end.
          pFstMatrix.alt.aor[i,j] = pFstMatrix.alt.aor[i,j] + Fst.ij
          
          # Add to pVector to calculate the median of all Fst values later.
          if (i != j) pVector = c(pVector, Fst.ij)
        }
        
        # Repeat the same process for varTable nad totTable. But we do not
        # calculate genome-wide Fst values by taking the average of ratios, but
        # we will take the ratio of averages (variance and pBar) at the end.
        if (!(is.null(varTable))){ 
          var.var.ij = var(c(P.var[i], P.var[j])) / 2
          varMatrix.var[i,j] = varMatrix.var[i,j] + var.var.ij
          
          pBar.var.ij = sum(W.i * P.var[i], W.j * P.var[j])
          pBarMatrix.var[i,j] = pBarMatrix.var[i,j] + pBar.var.ij
        }
        
        if (!(is.null(totTable))){ 
          var.tot.ij = var(c(P.tot[i], P.tot[j])) / 2
          varMatrix.tot[i,j] = varMatrix.tot[i,j] + var.tot.ij
          
          pBar.tot.ij = sum(W.i * P.tot[i], W.j * P.tot[j])
          pBarMatrix.tot[i,j] = pBarMatrix.tot[i,j] + pBar.tot.ij
        }
      }
      if (!(is.null(altTable))){
        # Take the median of the Fst values calculated for subpopulation i.
        pVector[i] = median(pVector[i:length(pVector)])
        pVector = pVector[1:i]
      }
    }
    # Store all pVectors in a matrix.
    if (!(is.null(altTable))){
      pFstMatrix.alt.median = rbind(pFstMatrix.alt.median, pVector) }
  }
  
  # Compute averages by dividing the sums by the the number of subpopulations.
  varMatrix.alt = varMatrix.alt / numRows
  varMatrix.var = varMatrix.var / numRows
  varMatrix.tot = varMatrix.tot / numRows
  
  pBarMatrix.alt = pBarMatrix.alt / numRows
  pBarMatrix.var = pBarMatrix.var / numRows
  pBarMatrix.tot = pBarMatrix.tot / numRows
  
  pFstMatrix.alt.aor = pFstMatrix.alt.aor / numRows
  
  # Computer genome-wide Fst values by taking the ratio of averages of 
  # variances and pBars. 
  pFstMatrix.alt.roa = varMatrix.alt / (pBarMatrix.alt * (1 - pBarMatrix.alt))
  pFstMatrix.var = varMatrix.var / (pBarMatrix.var * (1 - pBarMatrix.var))
  pFstMatrix.tot = varMatrix.tot / (pBarMatrix.tot * (1 - pBarMatrix.tot))
  
  # Assign appropriate row and column names for plotting purposes.
  colnames(varMatrix.alt) = rownames(varMatrix.alt) = populations
  colnames(varMatrix.var) = rownames(varMatrix.var) = populations
  colnames(varMatrix.tot) = rownames(varMatrix.tot) = populations
  
  colnames(pFstMatrix.alt.roa) = rownames(pFstMatrix.alt.roa) = populations
  colnames(pFstMatrix.alt.aor) = rownames(pFstMatrix.alt.aor) = populations
  colnames(pFstMatrix.var) = rownames(pFstMatrix.var) = populations
  colnames(pFstMatrix.tot) = rownames(pFstMatrix.tot) = populations
  
  if (!(is.null(altTable))){
    colnames(pFstMatrix.alt.median) = populations
    rownames(pFstMatrix.alt.median) = 
      paste(coordinates$CHR, ":", coordinates$POS, sep = "")
  }
  
  #___________________________ANALYZE SUPERPOPULATIONS____________________________
  if (!(is.null(altTable) & is.null(super))){ 
    cat("Analyzing super-populations...\n")
    # Create the superpopulations table from altTable.
    superNames = super$superpopulations
    super.s = super$super.s
    super.e = super$super.e
    sizes.super = c()
    for (pop in 1:length(superNames)){
      sizes.super = c(sizes.super, sum(sizes[super.s[pop]:super.e[pop]]))
    }

    altTable.super = matrix(ncol = 5, nrow = nrow(altTable))
    colnames(altTable.super) = superNames
    # Change NA's to 0's.
    altTable[is.na(altTable)] = 0
    for (i in 1:length(superNames)){
      altTable.super[, i] = apply(altTable[, super.s[i]:super.e[i]], 
                                  MARGIN = 1, FUN = sum) 
    }
    
    # Same procedure as subpopulations from here on...
    varMatrix.super = pBarMatrix.super = pFstMatrix.super =
      matrix(data=0, nrow=length(superNames), ncol=length(superNames))

    for (r in 1:numRows){
      row.super = as.numeric(altTable.super[r,])
      
      row.super[is.na(row.super)] = 0
      P.super = row.super / (sizes.super * 2)
      for (i in 1:length(superNames)){
        for (j in 1:length(superNames)){
          size = sizes.super[i] + sizes.super[j]
          W.i = sizes.super[i] / size
          W.j = sizes.super[j] / size
          # Compute variance.
          var.super.ij = var(c(P.super[i], P.super[j])) / 2
          varMatrix.super[i, j] = varMatrix.super[i, j] + var.super.ij
          # Compute pBar
          pBar.super.ij = sum(W.i * P.super[i], W.j * P.super[j])
          pBarMatrix.super[i, j] = pBarMatrix.super[i, j] + pBar.super.ij
        }
      }
    }

    varMatrix.super = varMatrix.super / numRows
    pBarMatrix.super = pBarMatrix.super / numRows
    
    pFstMatrix.super = 
      varMatrix.super / (pBarMatrix.super * (1 - pBarMatrix.super))
    colnames(pFstMatrix.super) = rownames(pFstMatrix.super) = superNames
  }
  #________________________________RETURN_______________________________________
  # Return a list of the matrices.
  results = list()
  results[[1]] = varMatrix.alt
  results[[2]] = varMatrix.var
  results[[3]] = varMatrix.tot
  
  results[[4]] = pBarMatrix.alt
  results[[5]] = pBarMatrix.var
  results[[6]] = pBarMatrix.tot
  
  results[[7]] = pFstMatrix.alt.roa
  results[[8]] = pFstMatrix.alt.aor
  results[[9]] = pFstMatrix.var
  results[[10]]= pFstMatrix.tot
  
  results[[11]] = pFstMatrix.alt.median
  
  results[[12]] = pFstMatrix.super
  names(results) = c("varMatrix.alt", "varMatrix.var", "varMatrix.tot",
                     "pBarMatrix.alt", "pBarMatrix.var", "pBarMatrix.tot",
                     "pFstMatrix.alt.roa", "pFstMatrix.alt.aor",
                     "pFstMatrix.var", "pFstMatrix.tot", 
                     "pFstMatrix.alt.median",
                     "pFstMatrix.super")
  return(results)
}
