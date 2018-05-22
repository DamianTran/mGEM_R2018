# Web scraping algorithm using XML parser and R url library

if(!require("XML")){
  install.packages("XML")
}
library(XML)

if(!require("RCurl")){ # We will be using RCurl for HTTPS security protocol support
  install.packages("RCurl")
}
library(RCurl)

# Example: getting synonyms to a word from thesaurus.com

findSynonyms = function(query){
  
  synonyms = readHTMLList(getURL(paste("http://www.thesaurus.com/browse/", query, "?s=t", sep = "")))
  # We paste our query into the variable field of the url
  # The URL, if protection is HTTP, can be pasted directly into readHTMLList()
  # However for flexibility, we wrap using getURL to ensure information can be extracted even from HTTPS domains
  
  # This will provide a garble of isolatable XML-type elements
  # We need to parse this list in order to extract word elements
  
  wordList = unique(synonyms[[6]]) # Printing the synonyms list to the console gives us the indication of where our words are located
  
  # To remove html header elements, we parse to remove any undue punctuation
  
  for(i in 1:length(wordList)){
    if(wordList[i] == query) wordList[i] = NA
    else if(nchar(wordList[i]) < 1) wordList[i] = NA
    else if((length(grep(",", wordList[i])) > 0) || (length(grep(";", wordList[i])) > 0) ||
            (length(grep("&", wordList[i])) > 0)){
      wordList[i] = NA;
    }
  }
  
  return(wordList[!is.na(wordList)]) # remove all of our values at the same time instead of reallocating the buffer upon deletion

}

pubMedSearch = function(query){
  
  # The information in the pubMed HTML is more nested, so we will get into the lower level
  # XML parser API
  
  tryCatch(
    {
      doc = htmlParse(getURL(paste("https://www.ncbi.nlm.nih.gov/pubmed/?term=", paste(query, collapse = "+"), sep = "")))
    },
    catch = function(){ 
      return(character(0)) 
    }
  )
  # We create an external pointer called "doc" to a structure that we create temporarily on the hard drive using htmlParse()
  termList = xpathSApply(doc, "//p", xmlValue) # This uses the XPath notation to acquire XML values hidden in the HTML node structure
  
  # Parse through the terms to remove any extraneous values
  for(i in 1:length(termList)){
    if(termList[i] == query) termList[i] = NA
    else if(termList[i] == "Similar articles ") termList[i] = NA
    else if(nchar(termList[i]) < 1) termList[i] = NA
    else if((length(grep("\\n", termList[i])) > 0) || (length(grep(";", termList[i])) > 0) ||
            (length(grep("&", termList[i])) > 0) || (length(grep("\\t", termList[i])) > 0)){
      termList[i] = NA;
    }
  }
  
  return(termList[!is.na(termList)])
  
}

# Now, we can chain this together with our geneGOSearch protocol to scrape google scholar and confirm if the gene
# is known to aggregate

getLiteratureAnnotations = function(genes, query = NA, collapse = FALSE){
  output = list()
  for(i in 1:length(genes)){
    cat("Searching web for query: ", genes[i], "\n")
    tryCatch(
      {
        if(!is.na(query) && length(query > 0)){
          terms = pubMedSearch(genes[i])
          if(collapse){ # AND search
            for(j in 1:length(query)){
              terms = terms[grep(query[j], terms)]
            }
          }
          else{ # OR search
            matchIdx = numeric(0)
            for(j in 1:length(query)){
              matchIdx = c(matchIdx, grep(query[j], terms)) 
            }
            terms = terms[unique(matchIdx)]
          }
          if(length(terms > 0)) output[[genes[i]]] = terms # Only add to the list if we have hits for these genes
        }
        else{
          terms = pubMedSearch(genes[i])
          if(length(terms) > 0) output[[genes[i]]] = pubMedSearch(genes[i])
        }
      },
      catch = function(){  } # Simply move on in the case of an error
    )
  }
  return(output)
}
