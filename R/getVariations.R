#script.R
getVariations <- function(populations, start=1, end=22){
  mainPath = getwd()
  gffPath = paste(mainPath, "/miRNA.gff3", sep = "")
  len = length(populations)
  
  # Read miRNA coordinates file.
  gff <- read.table(gffPath, sep="\t", quote="")
  gff <- gff[gff[,3]=="miRNA_primary_transcript", c(1,4,5,9)]
  colnames(gff) <- c("CHR", "START", "END", "INFO")

#_______________________________________________________________________________
  findVariations <- function(population, start, end){
    vcfPath <- paste(mainPath, "/VCF/", population, sep = "")
    resultsFile <- paste(vcfPath, "/results.tsv", sep = "")
    
    results <- data.frame(CHR = character(0), NAME = character(0), 
                          POS = character(0), ALT = character(0), 
                          altFREQ = character(0), varFREQ = character(0),
                          totFREQ = numeric(0), SIZE = numeric(0), 
                          stringsAsFactors = FALSE)
    if (start == 1){
      write.table(results, file = resultsFile, sep = "\t", quote = FALSE)
    }
    
    for (chr in start:end){
  
      results = results[0,]
      
      cat("Analyzing ", pop, " chromosome ", chr, "... ", sep = "")
      vcf <- read.table(paste(vcfPath, "/chr", toString(chr), ".vcf.gz", sep = ""))
      cat("VCF File Loaded.\n")
      
      mirna <- gff[gff$CHR==paste("chr", chr, sep = ""),]
      if (nrow(mirna) != 0){
        for (i in 1:nrow(mirna)){
          mirStart <- mirna$START[i]
          mirEnd <- mirna$END[i]
          v <- vcf[(vcf[,2]>=mirStart & vcf[,2]<=mirEnd),]
          #print(v)
          if (nrow(v) > 0){
            mirName<-strsplit(toString(mirna$INFO[i]), ";")[[1]][3]
            name <- sub("Name=", "", mirName)
            for (r in 1:nrow(v)){
              totFreq <- 0
              varFreq <- c(NA)
              altFreq <- c(NA)
              for(c in 10:ncol(v)){
                if (chr %in% c(1:22, "X") & v[r,c] != "0|0"){
                  one <- as.numeric(strsplit(toString(v[r,c]), "|")[[1]][1])
                  two <- as.numeric(strsplit(toString(v[r,c]), "|")[[1]][3])
                  
                  totFreq <- totFreq + 1
                  if(one != 0){
                    # Add to varFreq & altFreq
                    if (is.na(varFreq[one])){
                      varFreq[one] <- 1
                      altFreq[one] <- 1
                    } else {
                      varFreq[one] <- varFreq[one] + 1
                      altFreq[one] <- altFreq[one] + 1
                    }
                  } 
                  
                  if (two != 0){
                    # Add to varFreq only if _two_ is a different allele than _one_.
                    if (two != one){
                      if (is.na(varFreq[two])){
                        varFreq[two] <- 1
                        altFreq[two] <- 1
                      } else {
                        varFreq[two] <- varFreq[two] + 1
                        altFreq[two] <- altFreq[two] + 1
                      }
                    } else {
                      altFreq[two] <- altFreq[two] + 1
                    }
                  }
                } 
              } 
              size <- ncol(v)-9
              variation <- totFreq/size * 100
              if (variation != 0){ 
                cat("\t", name, ":\n", "\t\t", variation, 
                    "% of genomes have a variation at chromosome ",
                    v[r,1], " position ", v[r,2], ".\n", sep = "")
                ref <- toString(v[r,4])
                alt <- strsplit(toString(v[r,5]), ",")[[1]]
                freq <- c()
                for (t in 1:length(varFreq)){
                  if (is.na(varFreq[t])){varFreq[t] = 0}
                  freq <- c(freq, varFreq[t] / size * 100)
                  cat("\t\t", freq[t], "% have a ", ref, "->", 
                      alt[t], " variation.\n", sep = "")
                }
                for (k in 1:length(alt)){
                  results[nrow(results) + 1,] <- c(v[r,1], name, v[r,2], 
                                                  toString(alt[k]), as.numeric(altFreq[k]),
                                                  as.numeric(varFreq[k]), totFreq, size)
                }
              }
            }
          }
        }
      }
      write.table(results, file = resultsFile, append=TRUE, sep = "\t", 
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
  }
#_______________________________________________________________________________
  
  for (pop in populations){
    findVariations(pop, start, end)
  }
}
