# Core package file that will act as a header to source all other files

if(!require("rstudioapi")){ # R Studio API gives access to some more native functions of the interface
  install.packages("rstudioapi")
}
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

source("GO.R")
source("UniProt.R")
source("webscrape.R")

# The FINAL FORM - this function will query genes below a certain size that match the terms/function you search for, and then
# scrape PubMed to find the literature annotations for them and see if the data matches reality

findVerifiedGenes = function(query, size, ...){
  genes = findGenes(query, size, ...)
  return(getLiteratureAnnotations(genes, query))
}