# FIRST - this package requires BLAST+ to be installed
# Link can be found using ?blast
# Download the executable for your OS
# Install it, and then RESTART RStudio
# Check the proper installation with Sys.which("blastn")

# Check and install packages

if(!require("devtools")){
  install.packages("devtools") # We require devtools for the install_github() function
}
library(devtools)

if(!require("rBLAST")){
  install_github("mhahsler/rBLAST") # Install link grabbed from the rBLAST github readme.md
}
library(rBLAST)

if(!require("seqinr")){ # Contains simple functions for reading and writing FASTA files
  install.packages("seqinr")
}
library(seqinr)

if(!require("RFLPtools")){ # Contains read.blast() function that is a customized read.delim function for BLAST output types 6 or 10
  install.packages("RFLPtools")
}
library(RFLPtools)

# Example
# seq = *nucleotide sequence*
# write.fasta(seq, "Test Sequence", "test.fa")
# system("blastn -db nr -query test.fa -remote -outfmt 6 -out test.txt")

# Function to send aggregated count sequences for mapping to nucleotide BLAST
# Warning: takes a long time

BLAST_counts = function(counts, filebase, taxonID = 9606){
  
  # Counts: aggregated counts produced by the aggregate() function
  # filebase: the base of the filename, which will have .fa and .txt extensions appended to it by this function
  # taxonID: the id of the reference species - all others will be filtered out
  
  inFILE = paste(filebase, ".fa", sep = "")
  outFILE = paste(filebase, ".txt", sep = "")
  
  write.fasta(counts[1][1], "SEQUENCE_1", file.out = inFILE)
  for(i in 2:length(counts$uniqueSeq)){
    write.fasta(counts$uniqueSeq[i], paste("SEQUENCE_", i, sep = ""), file.out = inFILE, open = "a") # Append extra sequences to the same file
  }
  
  # All sequences are now appended, we will send a batch request to NCBI which will be processed faster than individual requests
  # due to the slow speed of network communication
  
  cat("Sending query with", nrow(counts), "sequences to NCBI Nucleotide BLAST...")
  system(paste("blastn -db nr -query ", inFILE, 
               " -remote -outfmt 6 -out ", outFILE, sep = ""))
  cat("\nResponse saved in", outFILE, '\n')
  
  # SUGGESTION: blast the most common sequences and be able to reverse-engineer the wild-type, 
  # then use the %identity to mathematically determine the %identity of other sequences 
  # based on their match to the wildtype - this will likely be much faster
  
  tryCatch({
    output = read.blast(outFILE)
    
    # Since mappings are returned in accession IDs of various nature, we use the MyGene servers
    # (From Umaseh's suggestion in MyGene.R) to discover the gene names
    # We will use the query function to determine which sequences map, and remove the ones that don't
    cat("Converting accession IDs\n")
    keep_index = numeric(0)
    for(i in 1:nrow(output)){
      tryCatch({
        check = query(output$subject.id[i])
        check$hits = check$hits[check$hits$taxid == taxonID,]
        if(nrow(check$hits) > 0){
          keep_index = c(keep_index, i)
          output$subject.id[i] = check$hits[1,]$name
        }
      },error = function(cond){ })
    }
    output = output[keep_index,] 
    write.table(output, sep = '\t', file = outFILE, col.names = FALSE, row.names = FALSE,
                quote = FALSE) # Replace the original with mapped conversions
    cat("SUCCESS - table saved to ", outFILE, '\n')
    return(output)
  }, error = function(cond){ 
    message(cond)
    return(NA)
  })
  # The returned dataframe can be used to process the input reads by match identity
}

# Function for filtering aggregated counts based on BLAST mappings to the expected sequence

BLAST_filter_counts = function(counts, name, BLASToutput, identityThreshold = 70){
  
  # counts: output produced by using aggregate() on a sequence matrix
  # name: the name of the gene of interest that we're looking for in this dataset
  # BLASToutput: an output produced by BLAST_counts on the input counts object
  # identityThreshold: how much must the query match the target gene (in %)
  
  BLASToutput = BLASToutput[grep(name, subject.id, ignore.case = TRUE),]
  mapNames = BLASToutput$query.id
  mapIdx = numeric(0)
  for(name in mapNames){
    mapIdx = c(mapIdx, as.numeric(substr(name, nchar(name), nchar(name))))
  }
  mapIdx = unique(mapIdx)
  counts = counts[mapIdx,]
}

# ALSO consider downloading the BLAST nucleotide database to your computer and process tasks
# locally, this is the way that frequent batch searches are typically done

