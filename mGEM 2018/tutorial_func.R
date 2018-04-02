# Sets the working directory - change this to your desired project directory
setwd("C:/Users/Damian/Documents/Tutorials/R/mGEM 2018/")

NGSgen = function(ncol, nrow, seqLength){
  randomMatrix = matrix(nrow = nrow, ncol = ncol)
  for(y in 1:nrow){
    for(x in 1:ncol){
      randomMatrix[y,x] = randSeq(seqLength)
    }
  }
  return(randomMatrix)
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