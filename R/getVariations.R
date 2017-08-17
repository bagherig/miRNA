#script.R
#' Finds the variations within the miRNA coordinates and writes the results
#' to a file named default 'results.tsv'.
#'
#' @param populations A vector of population names.
#' @param start The start chromosome number.
#' @param end The end chromosome number.
#'
#' @examples getVariations(c("ACB", "BEB"), 1, 5)
#' 
getVariations <- function(populations, start=1, end=22, filename="results"){
  # Error checking:
  if (!(start %in% 1:22 & end %in% 1:22)){ 
    stop("'start' and 'end' must be a value between 1 to 22.")
  }
  if (start > end){
    stop("'start' cannot be larger than 'end'")
  }
  
  mainPath = getwd()
  gffPath = file.path(mainPath, "miRNA.gff3")
  len = length(populations)
  
  # Read miRNA coordinates file.
  gff = read.table(gffPath, sep="\t", quote="")
  gff = gff[gff[, 3] == "miRNA_primary_transcript", c(1,4,5,9)]
  colnames(gff) = c("CHR", "START", "END", "INFO")

#_______________________________________________________________________________
  findVariations <- function(population){
    vcfPath = file.path(mainPath, "VCF", population)
    filename = paste(filename, ".tsv", sep = "")
    resultsFile = file.path(vcfPath, filename)
    
    # Data frame to store results.
    results = data.frame(CHR = character(0), NAME = character(0), 
                         POS = character(0), ALT = character(0), 
                          altFREQ = character(0), varFREQ = character(0),
                          totFREQ = numeric(0), SIZE = numeric(0), 
                          stringsAsFactors = FALSE)
    # If we are staring from chr 1, overwrite the file and write column names.
    if (start == 1){
      write.table(results, file = resultsFile, sep = "\t", quote = FALSE)
    }
    
    for (chr in start:end){
      # Clear results. 
      # We will write the results to the file chromosome by chromosome.
      results = results[0,]
      
      cat("Analyzing ", pop, " chromosome ", chr, "... ", sep = "")
      # Read this chromosome's VCF file.
      vcf = read.table(paste(vcfPath, "/chr", toString(chr), 
                             ".vcf.gz", sep = ""))
      cat("VCF File Loaded.\n")
      
      # Subset miRNA coordinates of this chromosome.
      mirna = gff[gff$CHR == paste("chr", chr, sep = ""), ]
      if (nrow(mirna) != 0){
        for (i in 1:nrow(mirna)){
          mirStart = mirna$START[i]
          mirEnd = mirna$END[i]
          # Subset variations within the miRNA coordinates range.
          vcf_filter = vcf[(vcf[, 2] >= mirStart & vcf[, 2] <= mirEnd), ]
          if (nrow(vcf_filter) > 0){
            # Extract miRNA name.
            mirName = strsplit(toString(mirna$INFO[i]), ";")[[1]][3]
            name = sub("Name=", "", mirName)
            for (r in 1:nrow(vcf_filter)){
              totFreq = 0      # total-variation frequency.
              varFreq = c(NA)  # variation frequency.
              altFreq = c(NA)  # SNP frequency.
              # Skip the first 9 columns of VCF file.
              for(c in 10:ncol(vcf_filter)){
                if (vcf_filter[r, c] != "0|0"){
                  # Seperate the two alleles.
                  one = as.numeric(strsplit(toString(vcf_filter[r, c]), 
                                            "|")[[1]][1])
                  two = as.numeric(strsplit(toString(vcf_filter[r, c]), 
                                            "|")[[1]][3])
                  totFreq = totFreq + 1
                  
                  if(one != 0){
                    # Add to varFreq & altFreq
                    if (is.na(varFreq[one])){
                      varFreq[one] = 1
                      altFreq[one] = 1
                    } else {
                      varFreq[one] = varFreq[one] + 1
                      altFreq[one] = altFreq[one] + 1
                    }
                  } 
                  
                  if (two != 0){
                    # Add to varFreq only if 'two' is a different allele 
                    # than 'one'.
                    if (two != one){
                      if (is.na(varFreq[two])){
                        varFreq[two] = 1
                        altFreq[two] = 1
                      } else {
                        varFreq[two] = varFreq[two] + 1
                        altFreq[two] = altFreq[two] + 1
                      }
                    } else {
                      altFreq[two] = altFreq[two] + 1
                    }
                  }
                } 
              } 
              
              size = ncol(vcf_filter) - 9 # Number of samples in the population.
              variation = totFreq / size * 100
              if (variation != 0){ 
                cat("\t", name, ":\n", "\t\t", variation, 
                    "% of genomes have a variation at chromosome ",
                    vcf_filter[r, 1], " position ", 
                    vcf_filter[r, 2], ".\n", sep = "")
                # Reference allele:
                ref = toString(vcf_filter[r, 4])  
                # Alternative alleles:
                alt = strsplit(toString(vcf_filter[r, 5]), ",")[[1]] 
                freq = c()  # Vector to store frequency of alleles. 
                for (t in 1:length(varFreq)){
                  if (is.na(varFreq[t])) { varFreq[t] = 0 }
                  freq = c(freq, varFreq[t] / size * 100)
                  cat("\t\t", freq[t], "% have a ", ref, "->", 
                      alt[t], " variation.\n", sep = "")
                }
                
                # Store each alternative allele in a seperate row.
                for (k in 1:length(alt)){
                  results[nrow(results) + 1, ] = c(vcf_filter[r,1], name, 
                                                   vcf_filter[r,2], 
                                                   toString(alt[k]), 
                                                   as.numeric(altFreq[k]),
                                                   as.numeric(varFreq[k]), 
                                                   totFreq, size)
                }
              }
            }
          }
        }
      }
      # Append the results of this chromosome to the the file.
      write.table(results, file = resultsFile, append=TRUE, sep = "\t", 
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
  }
#_______________________________________________________________________________
  
  for (pop in populations){
    findVariations(population = pop)
  }
}
