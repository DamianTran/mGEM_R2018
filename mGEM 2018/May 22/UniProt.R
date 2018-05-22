# We begin by setting the working directory to the parent directory of this R script

if(!exists("KEY_VAR")){
  KEY_VAR = data.frame(cbind(character(0), character(0)), stringsAsFactors = FALSE)
  colnames(KEY_VAR) = c("dataset", "key")
}

# Import the UniProt database
if(!require("UniProt.ws")){
  source("https://bioconductor.org/biocLite.R") 
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

queryUniProt = function(keys, columns, keytype = NA, dataset = humanUP_DB){
  
  # Example keys: TP53, BDNF, AHR
  # Example columns: "ENTREZ_GENE", "SEQUENCE"
  # Dataset must be an S4 SQL access object like the UniProt.ws class
  
  datasetName = deparse(substitute(dataset))
  if(nchar(datasetName) < 1){
    cat("Invalid dataset used\n")
    return(data.frame())
  }
  
  cat("Using dataset: ", datasetName, "\n") 
  # It's good habit when automating processes to make the machine read out
  # As much information as possible in case something goes wrong

  # Since R uses a parser to run a script during run-time
  # Every typed variable can be deparsed to its equivalent string
  # We will use this to have our algorithm learn from the user, and memorize the proper key to access desired data
  
  if(length(grep(datasetName, KEY_VAR$dataset)) > 0){
    keytype = as.character(KEY_VAR[grep(datasetName, KEY_VAR$dataset),][[2]])
    # If the keytype has been learned, no need to search for it
  }
  
  if((length(keytype) < 1) || is.na(keytype)){
    keySearch = keytypes(dataset); # Get all the keytypes in the dataset
    #Since we don't know the type of key we need, we will write an automated trial-error algorithm to find it for us
    
    for(i in 1:length(keySearch)){
      # Exceptions are handled in R using the tryCatch function which adopts Java-like syntax (can add a finally clause, etc.)
      
      tryCatch( # Example of exception handling in R
        {
          cat("Checking mapping to key ", keySearch[i], "...\n")
          selection = select(dataset, keys, columns, keySearch[i])
          
          assign("KEY_VAR", rbind(KEY_VAR, cbind(datasetName, keySearch[i])), .GlobalEnv)
          # Every assignment in R is scoped, therefore we must specify the environment that we are assigning our variable within
          
          colnames(KEY_VAR) = c("dataset", "key")
          
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
    cat("Using key type: ", keytype, "\n")
    # The above is a time-costly algorithm, so once we find the proper mapping we can just save it or memorize it
    # and then use it for future function calls to retreive the correct column in one go
    return(select(humanUP_DB, keys, columns, keytype))
  }

}

# To find out what key types are available:
# keytypes(humanUP_DB)

# To find out what columns you can search for:
# columns(humanUP_DB)

# Finally, we culminate our efforts in a single convenient function

findGeneSubset = function(terms, columns = c("DATABASE(PFAM)", "SEQUENCE", "DOMAINS", "GO"), size = NA, keytype = NA, ...){
  output = queryUniProt(geneGOSearch(terms, ...), columns, keytype)
  if(!is.na(size)) return(output[nchar(output$SEQUENCE) < size,])
  return(output)
}

findGenes = function(terms, size = NA, keytype = NA, ...){
  output = queryUniProt(geneGOSearch(terms, ...), c("SEQUENCE"), keytype);
  if(!is.na(size)) return(output[nchar(output$SEQUENCE) < size,][[1]])
  return(output[[1]])
}

# Fastest implementation like so:
# output = findGenes(c("protein", "aggregat"), c("DATABASE(PFAM)", "SEQUENCE", "GO"), 500/3, keytype = "GENECARDS")




