if (!require("pacman")) install.packages("pacman")
pacman::p_load(mygene)


# Description: Returns data from MyGene.info R Client (Ensembl, PubMed, UniProt, NCBI refseq, etc.)
# Arguments:
#   genes: vector or list, or string of comma-separated gene symbols
# Example: getMyGeneData("BRCA1,BRCA2")
getMyGeneData = function(genes) {
  output = queryMany(genes, scopes="symbol", species="human", fields="all")
  
  # TODO: Potentially get more data by directly querying http://mygene.info/tryapi/ (MyGene REST API)
  return(output)
}
