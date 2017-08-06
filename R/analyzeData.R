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
analyzeData <- function(populations, sizes, tables){
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
  coordinates = as.data.frame(tables["coordinates"])
  names(coordinates) = c("CHR", "POS")
  altTable = varTable = totTable = NULL
  if ("altTable" %in% names(tables)){ 
    altTable = as.data.frame(tables["altTable"])}
  if ("varTable" %in% names(tables)){ 
    varTable = as.data.frame(tables["varTable"])}
  if ("totTable" %in% names(tables)){ 
    totTable = as.data.frame(tables["totTable"])}
  
  numRows = max(nrow(altTable), nrow(varTable), nrow(totTable))
  print(nrow(altTable))
  print(nrow(varTable))
  print(nrow(totTable))
  print(numRows)
  # Go through every row of the table.
  for (r in 1:numRows){
    # Print the row number.
    cat("Analyzing row", r, "\n")
    # Clear pVector for each row.
    pVector = c()
    
    if (!(is.null(altTable))){ 
      row.alt = as.numeric(altTable[r,]) # Read the current row.
      row.alt[is.na(row.alt)] = 0   # switch NA values to 0.
      P.alt = row.alt/(sizes*2)     # calculate the frequencies of this SNP.
    }
    if (!(is.null(varTable))){ 
      row.var = as.numeric(varTable[r,]) # Read the current row.
      row.var[is.na(row.var)] = 0   # switch NA values to 0.
      P.var = row.var/sizes         # calculate the frequencies of this SNP.
    }
    if (!(is.null(totTable))){ 
      row.tot = as.numeric(totTable[r,]) # Read the current row.
      row.tot[is.na(row.tot)] = 0   # switch NA values to 0.
      P.tot = row.tot/sizes         # calculate the frequencies of this SNP.
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
        #if (is.nan(Fst.ij)) Fst.ij = 0
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
  # Return a list of the matrices.
  results = list()
  results["varMatrix.alt"] = varMatrix.alt
  results["varMatrix.var"] = varMatrix.var
  results["varMatrix.tot"] = varMatrix.tot
  
  results["pBarMatrix.alt"] = pBarMatrix.alt
  results["pBarMatrix.var"] = pBarMatrix.var
  results["pBarMatrix.tot"] = pBarMatrix.tot
  
  results["pFstMatrix.alt.roa"] = pFstMatrix.alt.roa
  results["pFstMatrix.alt.aor"] = pFstMatrix.alt.aor
  results["pFstMatrix.var"] = pFstMatrix.var
  results["pFstMatrix.tot"] = pFstMatrix.tot

  results["pFstMatrix.alt.median"] = pFstMatrix.alt.median
  
  return(results)
}
