# Core package file that will act as a header to source all other files

if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstudioapi)

if(!require("writexl")){ # R Studio API gives access to some more native functions of the interface
  install.packages("writexl")
}

library(writexl)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

source("GO.R")
source("UniProt.R")
source("webscrape.R")
source("BioMart.R")
source("Entrez.R")
source("MyGene.R")

# The FINAL FORM - this function will query genes below a certain size that match the terms/function you search for, and then
# scrape PubMed to find the literature annotations for them and see if the data matches reality

findVerifiedGenes = function(query, size, ...){
  genes = findGenes(query, size, ...)
  # df <- data.frame(matrix(unlist(getLiteratureAnnotations(genes, query))))
  tmp <- write_xlsx(getLiteratureAnnotations(genes, query),path = "Search Results.xlsx")
}