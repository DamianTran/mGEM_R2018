# mGEM Dry Lab R Workshop - May 7 2018
if (!require("pacman")) install.packages("pacman")
pacman::p_load(RSQLite, DBI, AnnotationForge, AnnotationDbi, org.Hs.eg.db, hgu95av2.db, GO.db)

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

# For simplicity, the chained call is re-written as a complete function
geneGOSearch = function(terms, ...){ # The ellipses carry over parameters to contained functions
  return(getGOGenes(subsetGOSearch(terms, ...))) # The ellipses here will pass the default parameters if none are provided
}


