# mergeData.R
#
#________________________________Merge Data_______________________________________
#' Merge results.tsv files into a single table. The .tsv file must have the 
#' predefined template. 
#'
#' @param populations A vector of population names to be merged.
#' @param alt logical; whether to return a table of "altFREQ" columns.
#' @param var logical; whether to return a table of "varFREQ" columns.
#' @param tot logical; whether to return a table of "totFREQ" columns.
#'
#' @return A list of tables.
#'
#' @examples mergeData(c("ACB", "BEB"))
#' 
mergeData <- function(populations, alt = TRUE, var = TRUE, tot = TRUE){
  # Check if at least one of alt, var, or tot is TRUE.
  if (!(alt | var | tot)){
    stop("At least one of alt, var, or tot must be TRUE.")
  }
  # Define the path to VCF folder.
  path = file.path(getwd(), "VCF")
  # Define the name of the file to be read.
  fileName = "results.tsv"
  # Variable for storing the merged table.
  mer = NULL
  # Open and merge the results file of each subpopulation.
  for (pop in populations){
    cat ("Merging ", pop, "...\n", sep = "")
    # Read the results.tsv file of this population.
    table <- read.table(file.path(path, pop, fileName), sep = "\t", 
                        header = TRUE, stringsAsFactors = FALSE)
    # Exclude any row that the total number of variations is less than 5 (since
    # they are insignificant). Also exclude column "SIZE".
    table <- table[table$totFREQ >= 5, 
                   c("CHR", "POS", "ALT", "altFREQ", "varFREQ", "totFREQ")]
    
    # Merge the tables by columns "CHR", "POS", and "ALT". Suppress warnigns for 
    # duplicate column names.
    if (is.null(mer)){
      mer <- table
    } else {
      suppressWarnings(mer <- merge(mer, table, all = TRUE, 
                                    by = c("CHR", "POS", "ALT")))
    }
    
    # Remove duplicated rows.
    mer = mer[!(duplicated(mer)), ]
  }
  
  # Seperate "altFREQ", "varFREQ", and "totFREQ" columns into seperate tables.
  tables = vector(mode = "list", length = 4)
  if (alt){ 
    # Select columns associcated with SNP frequency.
    altTable = mer[, seq(4, ncol(mer), 3)]
    altTable[is.na(altTable)] = 0
    colnames(altTable) = populations
    tables[[1]] = altTable
  }
  
  if (var){ 
    # Select columns associcated with variation frequency.
    varTable = mer[, seq(5, ncol(mer), 3)] 
    varTable[is.na(varTable)] = 0
    colnames(varTable) = populations
    tables[[2]] = varTable
  }
  
  if (tot){ 
    # Select columns associcated with total-variation frequency.
    totTable = mer[, c(1, 2, seq(6, ncol(mer), 3))]
    # totTable can still contain duplicated rows. Remove them.
    totTable = totTable[!(duplicated(totTable)), 3:ncol(totTable)]
    totTable[is.na(totTable)] = 0
    colnames(totTable) = populations
    tables[[3]] = totTable
  }
  
  # Table of "CHR" and "POS"
  tables[[4]] = mer[, c("CHR", "POS", "ALT")]
  names(tables) = c("altTable", "varTable", "totTable", "coordinates")
  return(tables)
}
