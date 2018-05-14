# We begin by setting the working directory to the parent directory of this R script

if(!require("rstudioapi")){ # R Studio API gives access to some more native functions of the interface
  install.packages("rstudioapi")
}
library(rstudioapi)

# We get the context of this document, and set the working directory to the path field
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Now that the home directory is known, we can source other files within it safely
source("GO.R")

# Import the UniProt database
if(!require("UniProt.ws")){
  biocLite("UniProt.ws")
}
library(UniProt.ws)

if(!exists("humanUP_DB")){
  availableUniprotSpecies(pattern = "Homo sapiens")
  humanUP_DB = UniProt.ws(9606)
}

# the UniProt.ws object is an interface object that can be constructed with the ID of the taxon of interest
# from OOP languages, the () operator on an object denotes "creation" or construction
# In effect, this action creates a class into humanUP_DB with the current species field set to "Homo sapiens"

# The select() function will allow us to get data from the dataset knowing the keys, column, and type of key that we want
# We will write a function that can simplify this search process for us

queryUniProt = function(keys, columns, keytype = NA){
  
  # Example keys: TP53, BDNF, AHR
  # Example columns: "ENTREZ_GENE", "SEQUENCE"
  
  if(is.na(keytype)){
    keySearch = keytypes(humanUP_DB); # Get all the keytypes in the dataset
    #Since we don't know the type of key we need, we will write an automated trial-error algorithm to find it for us
    
    for(i in 1:length(keySearch)){
      
      # Exceptions are handled in R using the tryCatch function which adopts Java-like syntax (can add a finally clause, etc.)
      
      tryCatch( # Example of exception handling in R
        {
          cat("Checking mapping to key ", keySearch[i], "...\n")
          selection = select(humanUP_DB, keys, columns, keySearch[i])
          # This key maps to our input, we can return the selection
          return(selection)
        },
        error = function(cond){ 
          message(cond)
          cat('\n')
          # No data for this key, keep searching
        })
      
    }
  }
  else{
    # The above is a time-costly algorithm, so once we find the proper mapping we can just save it or memorize it
    # and then use it for future function calls to retreive the correct column in one go
    return(select(humanUP_DB, keys, columns, keytype))
  }

}

# To find out what key types are available:
# keytypes(humanUP_DB)

# To find out what columns you can search for:
# columns(humanUP_DB)

tutorial_may14 = function(){

  # Now, we can integrate our GO gene search with this function:
  geneSearch = queryUniProt(geneGOSearch(c("protein", "aggregat")), c("KEGG", "DATABASE(PFAM)", "SEQUENCE"), "GENECARDS")
  
  # Get the subset of this search with small proteins with a sequence length < 500 nucleotides (each protein is 3 nucleotides)
  smallProt = geneSearch[nchar(geneSearch$SEQUENCE) < 500/3, ]

}

# Finally, we culminate our efforts in a single convenient function

findGenes = function(terms, columns = c("DATABASE(PFAM)", "SEQUENCE"), size = NA, keytype = NA, ...){
  output = queryUniProt(geneGOSearch(terms, ...), columns, keytype)
  if(!is.na(size)) return(output[nchar(output$SEQUENCE) < size,])
  return(output)
}

# Fastest implementation like so:
# output = findGenes(c("protein", "aggregat"), c("DATABASE(PFAM)", "SEQUENCE", "GO"), 500/3, keytype = "GENECARDS")




