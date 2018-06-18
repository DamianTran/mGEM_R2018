# Demo for simplified NGS processing

NGSdemo = function(){
  
  if(nchar(Sys.which("blastn")) < 1){
    cat("Could not detect blastn installation - please install or reinstall blast+, or restart RStudio\n")
    return()
  }

  # Amyloid-beta cDNA sequence
  ABseq = paste("GATGCAGAATTCCGACATGACTCAGGAT",
          "ATGAAGTTCATCATCAAAAATTGGTGTT",
          "CTTTGCAGAAGATGTGGGTTCAAACAAA",
          "GGTGCAATCATTGGACTCATGGTGGGCGG",
            "TGTTGTCATAGCGTAACATCATC",
          "ACCATCACCACTAA", sep = "")
  
  cat("Generating simulated NGS dataset")
  assign("NGSdemoset", simNGSmatrix(10,10, ABseq, 0.1, 0.9), envir = globalenv()) # Assign to global env
  # If you have time, crank the number of reads up higher - right now we are at 1m reads
  cat("\nAggregating counts\n")
  frequencies = aggregate(NGSdemoset, 0.7, 1.3) # Aggregate and size-select
  frequencies = frequencies[1:1000,] # Take the top 1000 because BLASTing 1m sequences would be absurd
  BLAST_check = BLAST_counts(frequencies, "NGSdemo_BLAST")
  cat("Filtering by BLAST match to \"Amyloid\"\n")
  frequencies = BLAST_filter_counts(frequencies, "Amyloid", 0.7) # BLAST to determine if sequences match the query
  
  # Normalize counts to log2
  frequncies$counts = log2(frequencies$counts)
  
  barplot(frequencies$counts, names.arg = frequencies$uniqueSeq, cex.names = 0.5, las = 2,
          space = 0,
          main = "Mutation Frequencies", ylab = "Counts") # Make the plot
  
  assign("NGSdemo_freq", frequencies, envir = globalenv())
  cat("Complete\n")
}
