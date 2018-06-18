if(!require("rstudioapi")){ # Load the R studio API and get the document context to set the working directory
  install.packages("rstudioapi")
}
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

source("NGSgen.R")

# Sample:
reads = simNGSmatrix(100,100,randSeq(20),0.1,0.9) # Generate a simulated NGS set
freq = aggregate(reads) # The entirety of read counts forms our null distribution

# The 95% confidence interval (CI) is 1.96 standard errors away from the mean
# The 99% CI is 2.58 standard errors away from the mean
# Therefore anything above or below that error is significantly deviant from the null distribution

st_dev = sd(freq$counts)
st_err = sd(freq$counts)/sqrt(length(freq$counts)) 
avg = mean(freq$counts)

# We will determine the p-value of an estimate (a count) from the mean by standard error approximation
# Our test is one-tailed: we are looking for significant elevations of a count from the null distribution

mutCount = freq$counts[2] # We won't use the wild-type for this example
z_score = (mutCount - avg)/st_err
p_value = pnorm(z_score, lower.tail = FALSE) # Upper tail only - anything lower will return p = 1

if(p_value <= 0.05){
  cat("The sequence ", freq$uniqueSeq[1], " is significantly enriched in the population\n")
}else{
  cat("The sequence ", freq$uniqueSeq[1], " is not significantly enriched in the population\n")
}
cat("Dataset avg: ", avg, "\nDataset SD: ", st_dev,
    "\nCount: ", mutCount, "\nTest score: ", z_score, 
    "\nP value: ", p_value, '\n\n')

mutCount = freq$counts[200] # A value lower on the sorted list will be unlikely to be significantly enriched
z_score = (mutCount - avg)/st_err
p_value = pnorm(z_score, lower.tail = FALSE)
if(p_value <= 0.05){
  cat("The sequence ", freq$uniqueSeq[1], " is significantly enriched in the population\n")
}else{
  cat("The sequence ", freq$uniqueSeq[1], " is not significantly enriched in the population\n")
}
cat("Dataset avg: ", avg, "\nDataset SD: ", st_dev,
    "\nCount: ", mutCount, "\nTest score: ", z_score, 
    "\nP value: ", p_value, '\n\n')

mutCount = freq$counts[500] # A value lower on the sorted list will be unlikely to be significantly enriched
z_score = (mutCount - avg)/st_err
p_value = pnorm(z_score, lower.tail = FALSE)
if(p_value <= 0.05){
  cat("The sequence ", freq$uniqueSeq[1], " is significantly enriched in the population\n")
}else{
  cat("The sequence ", freq$uniqueSeq[1], " is not significantly enriched in the population\n")
}
cat("Dataset avg: ", avg, "\nDataset SD: ", st_dev,
    "\nCount: ", mutCount, "\nTest score: ", z_score, 
    "\nP value: ", p_value, '\n\n')

plotFreq(reads)