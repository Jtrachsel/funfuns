#' GO term enrichment from non-model organisms
#'
#' @param Int_genes Genes/proteins that are enriched in some way/group
#' @param mapping_file Go term mapping file containing all genes/proteins in your experiment
#' @param ont Which ontology should be used? CC MF or BP
#' @param algor which algorithm should be used in the test? see topGO package
#' @param statistic which test statistic should be used? see topGO package
#'
#' @return returns a table of the test results for each GO term detected in your int_genes table
#' @export
#'
#' @examples None yet.
topGO_NonModel <- function(Int_genes, mapping_file, ont='BP', algor = 'elim', statistic='Fisher'){

  require(topGO)

  coreGenes <- Int_genes

  geneID2GO <- readMappings(mapping_file)
  geneNames <- names(geneID2GO)

  # Get the list of genes of interest
  myInterestingGenes <- coreGenes$accno
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  head(geneList)

  GOdata <- new("topGOdata", ontology = ont, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  # Run topGO with elimination test
  resultTopGO.elim <- runTest(GOdata, algorithm = algor, statistic = statistic )
  allRes <- GenTable(GOdata, pval = resultTopGO.elim,
                     orderBy = "pval",
                     topNodes = length(GOdata@graph@nodes),
                     numChar=1000)
  allRes <- allRes %>% mutate(ont=ifelse(ont=='BP', 'Biological Process',
                                         ifelse(ont=='MF', 'Molecular Function', "Cellular Component"))) %>%
    mutate(GO_aspect = ont,
           algorithm = algor,
           statistic = statistic) %>% dplyr::select(-ont)
  return(allRes)
  #write.table(allRes, file = "Lact_vivo_topGO_BP_results.txt", sep = "\t", quote = F, col.names = T, row.names = F)

}
