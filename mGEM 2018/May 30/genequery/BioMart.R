if (!require("pacman")) install.packages("pacman")
pacman::p_load(mygene, biomaRt)

# Description: Returns data from BioMart R Client (GO, Ensembl, NCBI refseq, etc.)
# Arguments:
#   genes: vector or list, or string of comma-separated gene symbols
#   attributes: vector of attributes that one wants to retrieve (= the output of the query)
# Example: getBioMartData(c('BRCA1', 'BRCA2'), c('chromosome_name', 'band'))
getBioMartData = function(genes, attributes = c('entrezgene')) {
  mygene_data = queryMany(genes, scopes="symbol", species="human", fields="ensembl.gene")
  ensembl_ids = mygene_data@listData[["ensembl.gene"]]
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  output = getBM(attributes=attributes, filters='ensembl_gene_id', values=ensembl_ids, mart=ensembl)
  
  # TODO: Potentially get more data by directly querying https://rest.ensembl.org/ (Ensembl REST API)
  return(output)
}