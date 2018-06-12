# Functions for simulated NGS data generation
# The NGSgen function generates a simulated output from a targetted NGS run
# ie. a specific region of DNA was amplified and what we are looking for are mutations inside this region

# Real sequencing files often come in the form of a *.BAM file containing a matrix of all the raw reads
# found in the sample by the Illumina HiSeq
#   These files are often multiple Gb in size with tens/hundreds of millions of strings

# In real NGS datasets, we need to account for a few factors:
#   Sample size (ie. number of reads)
#   The size of the gene in base pairs
#   The background noise
#   Contamination and foreign sequences

# Generate a noise matrix
randSeqMatrix = function(ncol, nrow, seqLength){
    randomMatrix = matrix(nrow = nrow, ncol = ncol)
    for(y in 1:nrow){
      for(x in 1:ncol){
        randomMatrix[y,x] = randSeq(seqLength)
      }
    }
    return(randomMatrix)
}

# Generate a perfect target sequence matrix
seqMatrix = function(ncol, nrow, seqLength){
  sequence = randSeq(seqLength)
  return (matrix(sequence, ncol, nrow))
}

# Function for random permutation of sequences to introduce noise
permutateMatrix = function(dataset, frequency, minSize, maxSize){
  
  # Dataset: input matrix of sequences
  # Frequency: 0 - 1 percent chance of introducing noise at each index
  # minSize, maxSize: the minimum and maximum size possibilities for the random sequences
  
  sizes = seq(minSize, maxSize)
  seqPerm = numeric(length(sizes))
  maxIndex = numeric(1)

  for(y in 1:nrow(dataset)){
    for(x in 1:ncol(dataset)){
      if(runif(1) < frequency){
        seqPerm = runif(length(sizes))
        maxIndex = which.max(seqPerm)
        dataset[y,x] = randSeq(sizes[maxIndex])
      }
    }
  }
  
  return(dataset)
}

# Introduce random truncations (ie. if DNA fragments due to degradation or sample handling)
permutateTruncate = function(dataset, frequency){
  rand = numeric(1)
  for(y in 1:nrow(dataset)){
    for(x in 1:ncol(dataset)){
      rand = runif(1)
      if(rand < frequency/2){ # Truncate from end
        seqPerm = runif(nchar(dataset[y,x])-1)
        maxIndex = which.max(seqPerm)
        dataset[y,x] = substr(dataset[y,x], 1, maxIndex)
      }
      else if(rand > 1.0-frequency/2){ #Truncate from beginning
        seqPerm = runif(nchar(dataset[y,x])-1)
        maxIndex = which.max(seqPerm)
        dataset[y,x] = substr(dataset[y,x], nchar(dataset[y,x])-maxIndex, nchar(dataset[y,x]))
      }
    }
  }
  
  return(dataset)
}

# Generate a single random sequence
randSeq = function(length){
  if(is.na(length) || (length < 1)) return(character(0))
  randString = character(length)
  randNum = numeric(1)
  for(i in 1:length){
    randNum = runif(1)
    if((randNum >= 0) & (randNum < 0.25)) randString[i] = "A"
    else if((randNum >= 0.25) & (randNum < 0.5)) randString[i] = "T"
    else if((randNum >= 0.5) & (randNum < 0.75)) randString[i] = "C"
    else randString[i] = "G"
  }
  return(paste(randString, collapse = ""))
}

# Generate random mutatnts in a sequence dataset
randMutate = function(dataset, frequency){
  colNum = ncol(dataset)
  rowNum = nrow(dataset)
  randNum = numeric(1)
  mutThreshold = 1-frequency
  for(y in 1:rowNum){
    for(x in 1:colNum){
      for(i in 1:nchar(dataset[y,x])){
        randNum = runif(1)
        if(randNum > mutThreshold){
          randNum = runif(1)
          if((randNum >= 0) & (randNum < 0.25)) substr(dataset[y,x], i, i) = "A"
          else if((randNum >= 0.25) & (randNum < 0.5)) substr(dataset[y,x], i, i) = "T"
          else if((randNum >= 0.5) & (randNum < 0.75)) substr(dataset[y,x], i, i) = "C"
          else substr(dataset[y,x], i, i+1) = "G"
        }
      }
    }
  }
  return(dataset)
}

# Finally, an access function to wrap all these permutators together
NGSgen = function(ncol, nrow, sequence, mutantFreq, quality){
  
  # ncol, nrow: x- and y- dimension lengths
  # length: target sequence length
  # mutantFreq: frequency of mutants
  # quality: inverse likelihood of noise effects (0-1)
  
  dataset = matrix(sequence, nrow, ncol)
  dataset = randMutate(dataset, mutantFreq)
  dataset = permutateMatrix(dataset, 1.0-quality, nchar(sequence)/2, 2*nchar(sequence))
  dataset = permutateTruncate(dataset, 1.0-quality)
  return(dataset)
}

# We can define classes in R to standardize processes such as this, and to provide more information to the user
setClass("simNGSmatrix",
         representation(wild_type = "character",
                        matrix = "matrix"),
         prototype = list(wild_type = character(),
                          matrix = matrix(0,0,0)))
# As it is, the "wild_type" and "matrix" elements must be accessed with the "@" operator

# Where sim_matrix is a simNGSmatrix class instance and sim_matrix@matrix would return the matrix it contains
# To make this class more portable, we define a method for the S4 class subsetting operator x$name
setMethod("$", "simNGSmatrix",
          function(x, name){
      slot(x, name)
})
setMethod("[", "simNGSmatrix", # We can also overload the indexing operator to access the matrix more easily
          function(x,i,j,...,drop = TRUE){ # Drop = TRUE allows us to use empty indexing such as x[1,]
            
  return(x@matrix[i,j]) # Use the base @ operator to access the "matrix" slot of x
            
  # Returing the value this way ensures our indexing is READ-ONLY so that data can never be modified by accident
  # by the user (ie. x[1,1] = "ATCG" is not possible, but we can read as in a = x[1,1])
})
setMethod("as.vector", "simNGSmatrix",
          function(x){
  return(as.vector(x@matrix))            
})
setMethod("==", signature(e1="simNGSmatrix",e2="ANY"),
          function(e1, e2){
  return(e1@matrix == e2)             
})
# Other operators can be overloaded for a custom class by setting methods for them like above

# This class is constructed using the R version of the new() operator
# ex. x = new("simNGSmatrix", wild_type = randSeq(n), matrix=);
# We can simplify this constructor call with an R function
simNGSmatrix = function(ncol, nrow, sequence, mutFreq, quality){
  return(new("simNGSmatrix", wild_type = sequence, matrix = NGSgen(ncol, nrow, sequence, mutFreq, quality)))
}

# Such that we can call a simNGSmatrix class instance with simNGSmatrix(10,10,randSeq(10),0.1,0.9)
# The wild-type sequence will be saved such that we can reference back to it in the future

# Aggregate the counts together in identical sequences
aggregate = function(dataset, lengthThreshold = 5){
  
  # Dataset: input sequence matrix
  #   Note: we have overloaded all the necessary operators to allow our simNGSmatrix class to fit right in here
  #   as if it were a normal matrix class, therefore we can pass a simNGSmatrix as an argument to this function
  # lengthThreshold: establish a minimum sequence length to remove short reads
  
  uniqueSeq = unique(as.vector(dataset))
  counts = numeric(length(uniqueSeq))
  for(i in 1:length(uniqueSeq)){
    counts[i] = sum(dataset == uniqueSeq[i]);
  }
  output = data.frame(uniqueSeq, counts, stringsAsFactors = FALSE)
  output = output[nchar(output$uniqueSeq) >= lengthThreshold,]
  return(output[order(output$counts, decreasing = TRUE),])
}

# A convenience function to automatically aggregate and plot the frequencies
plotFreq = function(simNGS, numSeq = 100, skipWT = TRUE,
                    main = "Mutation Frequencies",
                    ylab = "Count"){
  
  # simNGS: a simulated NGS dataset class, or matrix of reads
  # numSeq: number of top sequences to display
  # skipWT: skip the wild-type sequence and display mutants on the plot only
  # main: plot title
  # ylab: y label of plot
  
  frequencies = aggregate(simNGS)[1:numSeq,]
  if(skipWT) frequencies = frequencies[frequencies$uniqueSeq != simNGS$wild_type,]
  plot.new() # Create a new plot
  par(oma=c(2*max(nchar(frequencies$uniqueSeq))/8,0,0,0)) # Set the plot window parameters
  barplot(frequencies$counts, names.arg = frequencies$uniqueSeq, cex.names = 0.5, las = 2,
          space = 0,
          main = main, ylab = ylab) # Make the plot
}

# Try this:
#   plotFreq(simNGSmatrix(100,100,randSeq(20),0.1,0.9))

# Some statistics:
# Now being able to aggregate the mutants, we need to be able to determine if a certain mutation is occurring beyond chance

# source("stats.R")

