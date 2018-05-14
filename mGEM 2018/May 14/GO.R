# mGEM Dry Lab R Workshop - May 7 2018

if(!require("RSQLite")){ # Check if library exists, if fails then install
  install.packages("RSQLite")
}
library(RSQLite) #import the library to the workspace

if(!require("DBI")){ # DBI is the main R interface to database servers
  install.packages("DBI")
}
library(DBI)

if(!require("AnnotationForge") || !require("AnnotationDbi") || !require("org.Hs.eg.db") || 
   !require("hgu95av2.db") || !require("GO.db")){ # pre-check to see if we need the BioConductor installer
  source("https://bioconductor.org/biocLite.R") # Source the bioconductor installer from their server
  biocLite() # Install the installer
}

if(!require("AnnotationForge")){
  biocLite("AnnotationForge") # High-level library containing a toolkit for annotation database creation and access
}
library(AnnotationForge)

if(!require("AnnotationDbi")){
  biocLite("AnnotationDbi") # AnnotationDbi allows for interface with ".db" extension annotation filetypes
}
library(AnnotationDbi)

if(!require("org.Hs.eg.db")){
  biocLite("org.Hs.eg.db") # Library containing ".db" file format information
}
library(org.Hs.eg.db)

if(!require("hgu95av2.db")){
  biocLite("hgu95av2.db") # database file containing gene symbols and GO terms (in addition to microarray probe mappings, not required)
}
library(hgu95av2.db)

if(!require("GO.db")){
  biocLite("GO.db") # GO term database containing terms and definitions
}
library(GO.db)

if(!exists("probe2GO")){ # This is essentially a check to see if this file has been previously sourced by the user
  cat(">> Datasets contained in human gene probe set: \n\n")
  ls("package:hgu95av2.db") # Show the dataframes contained within the hgu95av2 database file
  cat("\n>> Datasets contained in the GO database: \n\n")
  ls("package:GO.db")
}

# Note that the command calls were separated to prevent the loading of all packages in the event of a single missing package
# biocLite (c("AnnotationDbi", "org.Hs.eg.db", "GO.db" )) + other packages would also work

GOdb = as.data.frame(GOTERM) # List of all GO terms and definitions

# Now, let's get the indices that contain words of interest

subsetGOSearch = function(terms, category = "Definition", collapse = TRUE){ 
  
  # We'll use a function to find indices associated with a variety of terms and collapse them
  
  # Terms: the search terms
  # Category: the column of the GO database to search from
  # Collapse: do we find rows that match all terms together (TRUE) or that contain any one of the terms we've given? (FALSE)
  
  idx = numeric(0)
  
  if(collapse){
    for(i in 1:length(terms)){
      if(length(idx) < 1) idx = grep(terms[i], GOdb[[category]]) # First term takes baseline subset
      else{
        idx = idx[match(grep(terms[i], GOdb[[category]]), idx)]# All subsequent terms are taken from the subset of the first
        idx = idx[!is.na(idx)] # match produces NA values when a match fails, so remove these from each round
      }
    }
  }
  else{
    for(i in 1:length(terms)){
      idx = c(idx, grep(terms[i], GOdb[[category]]));
    }
    idx = unique(idx)
    idx = idx[order(idx)]
  }
  
  return(GOdb[idx,]) # Return the subset of the GO database matching the search query
}

GOsubset = subsetGOSearch(c("protein", "aggregat")) # Example search will give GO terms associated with proteins OR aggregation
GOSpecificSubset = subsetGOSearch(c("protein", "ribosome"), "Term", TRUE) # Second search will give GO terms associated with protein AND ribosome

# Now a second function to map the GO terms to human genes

probe2GO = as.data.frame(hgu95av2GO) # Mapping of probes to GO terms
probe2Symbol = as.data.frame(hgu95av2SYMBOL) # mapping of probes to official (Hugo) gene symbols

getGOGenes = function(GOset){ # Pass the entire go set
  
  # Function requires probe2GO and probe2Symbol
  # GOset: the set containing a "go_id" field
  
  matchIdx = match(GOset$go_id, probe2GO$go_id);
  matchIdx = unique(matchIdx[!is.na(matchIdx)])
  
  probeSubset = probe2GO[matchIdx,]
  
  matchIdx = match(probeSubset$probe_id, probe2Symbol$probe_id)
  matchIdx = matchIdx[!is.na(matchIdx)]
  
  output = probe2Symbol[matchIdx,]
  
  return(unique(output$symbol[!is.na(output$symbol)])) # Return a character array containing the gene names that match these GO terms

}

subsetGenes = getGOGenes(GOsubset)
specificGenes = getGOGenes(GOSpecificSubset)

# ta-da!  Future implementations should include a condensed table with the gene names side-by-side to the GO annotations
# for a more complete reference that can be analyzed and checked for validity

# Here, we chain the calls together, getting the genes from the result of the GO subset search
interestGenes = getGOGenes(subsetGOSearch(c("protein", "aggregat"), "Definition", TRUE))

#Analog

# For simplicity, the chained call is re-written as a complete function
geneGOSearch = function(genes, ...){ # The ellipses carry over parameters to contained functions
  return(getGOGenes(subsetGOSearch(genes, ...))) # The ellipses here will pass the default parameters if none are provided
}


