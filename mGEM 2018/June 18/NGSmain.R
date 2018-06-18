if(!require("rstudioapi")){ # Load the R studio API
  install.packages("rstudioapi")
}
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set the working directory

source("NGSgen.R")
source("NGSprocess.R")
source("BLAST.R")
source("MyGene.R")
source("NGSdemo.R")