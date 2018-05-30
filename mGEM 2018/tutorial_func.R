# Sets the working directory - change this to your desired project directory
setwd("C:/Users/Umaseh/Documents/GitHub/mGEM_R2018/mGEM 2018/")

NGSgen = function(ncol, nrow, seqLength = NA, seq = NA){
  if(is.na(seq) & !is.na(seqLength)){
    randomMatrix = matrix(nrow = nrow, ncol = ncol)
    for(y in 1:nrow){
      for(x in 1:ncol){
        randomMatrix[y,x] = randSeq(seqLength)
      }
    }
    return(randomMatrix)
  }
  return (matrix(seq, ncol, nrow))
}

randSeq = function(length){
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

aggregate = function(dataset){
  uniqueSeq = unique(as.vector(dataset))
  counts = numeric(length(uniqueSeq))
  for(i in 1:length(uniqueSeq)){
    counts[i] = sum(dataset == uniqueSeq[i]);
  }
  return(data.frame(uniqueSeq, counts))
}