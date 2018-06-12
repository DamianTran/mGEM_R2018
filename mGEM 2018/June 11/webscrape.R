# Web scraping algorithm using XML parser and R url library

if (!require("pacman")) install.packages("pacman")
pacman::p_load(XML, RCurl)

if(!require("stringr")){
  install.packages("stringr")
}
library(stringr)

instr <- function(str1,str2,startpos=1,n=1){
  aa=unlist(strsplit(substring(str1,startpos),str2))
  if(length(aa) < n+1 ) return(0);
  return(sum(nchar(aa[1:n])) + startpos+(n-1)*nchar(str2) )
}
# Example: getting synonyms to a word from thesaurus.com

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

findSynonyms = function(query){
  
  synonyms = readHTMLList(getURL(paste("http://www.thesaurus.com/browse/", query, sep = "")))
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
      doc = getURL(paste("https://www.ncbi.nlm.nih.gov/pubmed/?term=", paste(query, collapse = "+"), "&dispmax=200&report=medline", sep = ""))
    },
    catch = function(){ 
      return(character(0)) 
    }
  )
  
  lines = unlist(strsplit(doc, "\n"), use.names=FALSE) #get lines of MEDLINE data
  PMIDs = sub("PMID- ", "", lines[str_detect(lines, "PMID-")]) #get PMID from MEDLINE
  URLs = sub("PMID- ", "https://www.ncbi.nlm.nih.gov/pubmed/", lines[str_detect(lines, "PMID-")]) #construct URL from PMID
  #titles =  sub("TI  - ", "", lines[str_detect(lines, "TI  -")]) #get titles from MEDLINE
  titles = c()
  abstracts = c()
  
  for (j in 1:length(lines)){
    if (str_detect(lines[j],"TI  -")|str_detect(lines[j],"AB  -")){
      i = j+1
      text = lines[j]
      while (!str_detect(lines[i]," -"))
      {
        text = paste(text,trim(lines[i]))
        i = i+1
      }
      if (str_detect(lines[j],"TI  -")){
        titles = c(titles, sub("TI  - ", "", text))
      }
      else{
        while (length(titles)>length(abstracts)+1){
          abstracts = c(abstracts,"No abstracts found")
        }
        abstracts = c(abstracts, sub("AB  - ", "", text))
      }
    }
  }
  
  while (length(titles)>length(abstracts)){
    abstracts = c(abstracts,"No abstracts found")
    }
  
  info = data.frame(titles, PMIDs,URLs, abstracts)
  colnames(info)= c("Title", "PMID","URL", "Abstract")
  return(info)
  
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
          info = pubMedSearch(genes[i])
          titles = unlist(info["Title"],use.names=FALSE)
          PMID = unlist(info["PMID"],use.names=FALSE)
          URLs = unlist(info["URL"],use.names=FALSE)
          abstracts = unlist(info["Abstract"],use.names=FALSE)
          if(collapse){ # AND search
            for(j in 1:length(query)){
              titles = titles[grep(query[j], titles)]
              PMID = PMID[grep(query[j], titles)]
              abstracts = abstracts[grep(query[j], titles)]
            }
          }
          else{ # OR search
            matchIdx = numeric(0)
            for(j in 1:length(query)){
              matchIdx = c(matchIdx, grep(query[j], titles)) 
            }
            titles = titles[unique(matchIdx)]
            PMID = PMID[unique(matchIdx)]
            URLs = URLs[unique(matchIdx)]
            abstracts = abstracts[unique(matchIdx)]
          }
          if(length(titles > 0)) 
          {
            result = data.frame(titles,PMID,URLs,abstracts)
            colnames(result)= c("Title", "PMID", "URL", "Abstract")
            output[[genes[i]]]=result
          }# Only add to the list if we have hits for these genes
        }
        else{
          output[[genes[i]]]=info
        }
      },
      catch = function(){  } # Simply move on in the case of an error
    )
  }
  return(output)
}
