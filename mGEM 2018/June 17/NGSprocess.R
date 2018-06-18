# Aggregate the counts together in identical sequences
aggregate = function(dataset, minThres = 0.7, maxThres = 1.3, minLength = 3){
  
  # Dataset: input sequence matrix
  #   Note: we have overloaded all the necessary operators to allow our simNGSmatrix class to fit right in here
  #   as if it were a normal matrix class, therefore we can pass a simNGSmatrix as an argument to this function
  # minThres, maxThres: establish sequence length factor thresholds to exclude non-mutation fragments
  # minLength: minimum number of base pairs to be considered - removes single nucleotides etc.
  
  uniqueSeq = unique(as.vector(dataset))
  counts = numeric(length(uniqueSeq))
  for(i in 1:length(uniqueSeq)){
    counts[i] = sum(dataset == uniqueSeq[i]);
  }
  
  output = data.frame(uniqueSeq, counts, stringsAsFactors = FALSE)
  output = output[(nchar(output$uniqueSeq) >= minLength),]

  highL = nchar(output[1,1])*maxThres
  lowL = nchar(output[1,1])*minThres
  output = output[((nchar(output$uniqueSeq) >= lowL) && (nchar(output$uniqueSeq) <= highL)),]
  
  return(output[order(output$counts, decreasing = TRUE),])
}

