if (!require("pacman")) install.packages("pacman")
pacman::p_load(rentrez)

# Description: Returns data for multiple genes from Entrez R Client (NCBI REST API)
#
# Warning: The NCBI will ban IPs that don't use EUtils within their user guidelines. In particular:
# Don't send more than three request per second (rentrez enforces this limit)
# If you plan on sending a sequence of more than ~100 requests, do so outside of peak times for the US
# For large requests use the web history method
#
# Arguments:
#   genes: vector or list of gene symbols
#   query: list of search terms
# Example: getEntrezData(c('SKP1', 'VAMP1'), c('protein', 'aggregat'))
getEntrezData = function(genes, query) {
  output = list()
  for (gene in genes) {
    output[[gene]] = getEntrezDataSingleGene(gene, query)
  }
  
  # TODO: Potentially get more data by directly querying https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ (NCBI REST API)
  return(output)
}

# Description: Returns data for single gene from Entrez R Client (NCBI REST API)
# Arguments:
#   genes: string of gene symbol
#   query: list of search terms
# Example: getEntrezDataSingleGene('SKP1', c('protein', 'aggregat'))
getEntrezDataSingleGene = function(gene, query) {
  wildcard_query = paste0(paste(query, collapse='* AND '), '*')
  search_query = paste('Homo[ORGN]', gene, wildcard_query, sep=' AND ')
  output = entrez_search(db='pubmed', term=search_query)
  return(output)
}