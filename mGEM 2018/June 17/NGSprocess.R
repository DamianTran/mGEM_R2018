if(!require("textreuse")){ # We will use the align_local() function here to check mutation positions
  install.packages("textreuse")
  library(textreuse)
}

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

# Function to port the align_local to DNA sequence alignment
alignDNA = function(template, test){
  
  # template: the sequence to be compared against
  # test: the sequence that we want to check the alignment of against the template
  
  # First we need to split the strings by spaces for the align_local function to work
  template = paste(strsplit(template, split = "")[[1]], collapse = " ")
  test = paste(strsplit(test, split = "")[[1]], collapse = " ")
  
  alignment_test = align_local(template, test)
  
  # Now we concatenate the string back together and check for edits indicated by the '#' character
  alignment_test$a_edits = paste(strsplit(alignment_test$a_edits, split = " ")[[1]], collapse = "")
  alignment_test$b_edits = paste(strsplit(alignment_test$b_edits, split = " ")[[1]], collapse = "")
  
  return(alignment_test)
  
}

# Function for aggregation of mutants by mutation position
aggregate_pos = function(counts, wild_type = NA){
  if(nrow(counts) < 1) return(NA)
  if(is.na(wild_type)) wild_type = counts$uniqueSeq[1]
  
  # TO DO: compare each sequence in the counts to the wild-type reference
  # Tally the location of the edits, denoted by '#', and add an extra column
  # to the counts dataframe containing the location of the edit
  
}

# TO DO: write another function for the annotation of mutations by their type
# ie. point mutation, silent mutation (will require a function that translates the DNA to amino acids)
# or frame shift (indel)

