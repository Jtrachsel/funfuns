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
#' @examples #None yet.
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



#' Cumulative distance to targets in a distance matrix
#'
#' @param tmat MATRIX representation of distance matrix
#' @param pattern that defines your targets, passed to grep
#'
#' @return returns a named vector of cumulative distance to your specified targets (for every entry)
#' @export
#'
#' @examples #Soon
cum_dist_to_targets <- function(tmat, pattern){

  dists_to_targs <- tmat[grep(pattern, rownames(tmat)),]
  dist_to_targets <- colSums(dists_to_targs)
  return(dist_to_targets)
}

#' Minimum distance to targets for each other entry in a distance matrix
#'
#' @param tmat MATRIX representation of distance matrix
#' @param pattern that defines your targets, passed to grep
#'
#' @return returns a named vector of minimum distance to your specified targets (for every entry)
#' @export
#'
#' @examples #soon
min_dist_to_targets <- function(tmat, pattern){

  dists_to_targs <- tmat[grep(pattern, rownames(tmat)),]
  min_dists <- apply(dists_to_targs,2,min)

  return(min_dists)
}

#' iteratively remove objects from a distance matrix and plot at each iteration
#'
#' @param in_dist input distance matrix, must be type dist
#' @param trymax passed to metaMDS function
#' @param iterations number of removal iterations to try
#' @param maxit passed to metaMDS function
#' @param autotransform passed to metaMDS function
#' @param exclusion_prob passed to quantile(), all observations with cumulative distances to targets (or total) above this quantile will be removed each iteration
#' @param parallel passed to metaMDS function
#' @param targ_pat passed to grep, defines the targets you are interested in (ones you dont want to remove).  Only used
#' for rem_type = 'cum_targ' and 'min'
#' @param rem_type one of 'total', 'cum_targ', and 'min'.  Specifies how to calculate distances for removal. total will pass cumulative distances to quantile, meaning the
#' objects with the greatest cumulative distance from all other objects will be removed first. 'cum_targ' considers cumulative distances to objects selected by your
#' targ_pat parameter, objects with greatest cumulative distances to your targets will be removed first. min considers the minimum distance to your targets, that is, the distance
#' from each object to the closest target. probably a good choice if your targets are similar to each other and you want to keep
#' all objects that are very similar to any of your targets.
#' @param nmds
#'
#' @return returns a big ol list
#' @export
#'
#' @examples #soon
iterative_dist_remove <- function(in_dist, trymax=50,
                     iterations = 20,
                     maxit=1000,
                     autotransform = FALSE,
                     exclusion_prob = 0.95,
                     parallel=8,
                     targ_pat = 'SX',
                     rem_type = 'total',
                     nmds=FALSE){
  require(ggrepel)
  require(vegan)
  require(Rtsne)
  require(outliers)
  final_results <- list()
  for (iter in 1:iterations){
    print(iter)

    # finds outliers based on pwdists
    if (iter ==1){
      dist <- in_dist
    } else {
      print(paste('using clean_mat from iter: ', iter -1))
      dist <- final_results[[(iter-1)]][[3]]
    }

    mainmat <- as.matrix(dist)
    # print(nrow(mainmat))
    # print(ncol(mainmat))
    if (rem_type == 'min'){
      dists_to_mine <- min_dist_to_targets(tmat = mainmat, pattern=targ_pat)
    }

    if (rem_type == 'total'){
      dists_to_mine <- rowSums(mainmat)
    }

    if (rem_type == 'cum_targ'){
      dists_to_mine <- cum_dist_to_targets(tmat = mainmat, pattern = targ_pat)

    }

    # print(hist(dists_to_mine))
    # using outliers is wrong because it will remove outliers with unusually small distances
    # as well as large, meaning very similar genomes are likely to be removed as well
    # using quantile like below will only remove genomes from the upper end of distances
    # badones <- outliers::scores(x = dists_to_mine, type = 't', prob = 0.95)
    # print(hist(dists_to_mine[badones], add=TRUE, color='blue'))
    # browser()
    badones <- dists_to_mine > quantile(dists_to_mine, probs = exclusion_prob)



    # badones <- outliers::scores(x = rowSums(mainmat), type = outlier_type, prob = outlier_prob)
    # mainmat <- as.matrix(test_dist)
    remove_me <- names(badones)[badones]
    rem_cols <- which(colnames(mainmat) %in% remove_me)
    rem_rows <- which(rownames(mainmat) %in% remove_me )
    mainmat <- mainmat[-rem_rows, -rem_cols]
    clean_dist <- as.dist(mainmat)

    # print(nrow(mainmat))
    # print(ncol(mainmat))

    if (nrow(mainmat) ==0){
      print(paste('None/All your data was removed at iteration, but maybe something else is happening', iter))
      #return(final_results)
      return(list(dists_to_mine, badones))
    } else{
      iter_results <- list()
      print(paste('number of genomes:', length(mainmat[1,])))

      ### tsne calc ###

      rtsne_test <- adapt_tsne(dist = dist)

      # rtsne_test <- Rtsne(dist, is_distance = TRUE, perplexity = 30, max_iter = 2000)
      tsne_points <- as.data.frame(rtsne_test$Y)
      tsne_points$genome <- attributes(dist)$Labels


      if (nmds == TRUE){
        ### NMDS calc
        NMDS <- metaMDS(dist, trymax = trymax, autotransform = autotransform, k=2, parallel = parallel, maxit=maxit)
        nmds <- as.data.frame(NMDS$points)
        nmds$genome <- rownames(nmds)
        all_points <- merge(nmds, tsne_points, by = 'genome')
        bads <- nmds[nmds$genome %in% remove_me,]
        ggplot(nmds, aes(x=MDS1, y=MDS2, label=genome, color=serovar)) + geom_point()
        p1 <- ggplot(nmds, aes(x=MDS1, y=MDS2)) + geom_point(data=bads, aes(x=MDS1, y=MDS2), color='purple') +
          geom_point() + ggtitle('NMDS: purple points will be removed...') + theme(legend.position = 'none')
        iter_results[[2]] <- p1
      } else{
        all_points <-tsne_points
      }

      bads <- all_points[all_points$genome %in% remove_me,]

      p2 <- all_points %>% ggplot(aes(x=V1, y=V2)) +
        geom_point() + theme(legend.position = 'none') + geom_point(data=bads, aes(x=V1, y=V2), color='purple') +
        ggtitle('TSNE: purple points will be removed...')

      # bads <- nmds[nmds$genome %in% remove_me,]
      # ggplot(nmds, aes(x=MDS1, y=MDS2, label=genome, color=serovar)) + geom_point()
      # p1 <- ggplot(nmds, aes(x=MDS1, y=MDS2)) + geom_point(data=bads, aes(x=MDS1, y=MDS2), color='purple') +
      # geom_point() + ggtitle('NMDS: purple points will be removed...') + theme(legend.position = 'none')

      #iter_results[[1]] <- nmds
      iter_results[[1]] <- all_points
      # iter_results[[2]] <- p1
      iter_results[[3]] <- clean_dist
      iter_results[[4]] <- p2
      # print(p1)
      print(p2)

      final_results[[iter]] <- iter_results
    }

  }
  return(final_results)
}


#' plots the output of iterative_NMDS
#'
#' @param iter_res the full list output by iterative_NMDS
#' @param iter the iteration you want to plot, should be an integer
#' @param pattern the pattern that will identify objects you are interested in, passed to grep
#' @param labels should we label things in the plots?
#' @param mysize what size do you want the points you are interested in to have?
#'
#' @return
#' @export
#'
#' @examples #soon
iter_plots <- function(iter_res, iter, pattern, labels = TRUE, mysize = 1){

  res <- iter_res[[iter]][[1]]
  mine <- res[grep(pattern = pattern, res$genome),]
  num_match <- length(mine$genome)
  tot <- length(res$genome)

  if (labels == TRUE){
    if ('MDS1' %in% colnames(res)){
      p.nmds <- res %>% ggplot(aes(x=MDS1, y=MDS2)) + geom_point(alpha=.7) +
        geom_point(data=mine, aes(x=MDS1, y=MDS2), fill='red', size = mysize, shape=21) +
        geom_text(data=mine, aes(x=MDS1, y=MDS2, label=genome)) +
        ggtitle("NMDS", subtitle = paste('total genomes:', tot, '\ngenomes matching pattern:', num_match, sep = ' '))
    } else {p.nmds <- NULL}
    p.tsne <- res %>% ggplot(aes(x=V1, y=V2)) + geom_point(alpha=.7) +
      geom_point(data=mine, aes(x=V1, y=V2), fill='red', size = mysize, shape=21) +
      geom_text(data=mine, aes(x=V1, y=V2, label=genome), alpha=.75) +
      ggtitle("TSNE", subtitle = paste('total genomes:', tot, '\ngenomes matching pattern:', num_match, sep = ' '))

  } else {
    if ('MDS1' %in% colnames(res)){
      p.nmds <- res %>% ggplot(aes(x=MDS1, y=MDS2)) + geom_point(alpha=.7) +
        geom_point(data=mine, aes(x=MDS1, y=MDS2), fill='red', size = mysize, shape=21) +
        ggtitle("NMDS", subtitle = paste('total genomes:', tot, '\ngenomes matching pattern:', num_match, sep = ' '))
    }else {p.nmds <- NULL}

    p.tsne <- res %>% ggplot(aes(x=V1, y=V2)) + geom_point(alpha=.7) +
      geom_point(data=mine, aes(x=V1, y=V2), fill='red', size = mysize, shape=21) +
      ggtitle("TSNE", subtitle = paste('total genomes:', tot, '\ngenomes matching pattern:', num_match, sep = ' '))

  }

  return(list(p.nmds, p.tsne))

}


#' wrapper for rtsne function that will set perplexity to 1 if an error occurs
#'
#' @param dist input distance matrix, must be type dist
#'
#' @return rtsne result
#' @export
#'
#' @examples #soon
adapt_tsne <- function(dist) {
  out <- tryCatch(
    {

      message("running Rtsne function...")
      Rtsne(dist, is_distance = TRUE, perplexity = 30, max_iter = 2000)


    },
    error=function(cond) {
      message('error trying to run Rtsne...')
      message("Here's the original error message:\n")
      message(cond)
      message('setting perplexity to 1')
      # Choose a return value in case of error
      Rtsne(dist, is_distance = TRUE, perplexity = 1, max_iter = 2000)
      # return(NA)
    },
    warning=function(cond) {
      message('warning trying to run Rtsne...')
      message("Here's the original warning message:\n")
      message(cond)
      message('setting perplexity to 1')
      # Choose a return value in case of warning
      Rtsne(dist, is_distance = TRUE, perplexity = 1, max_iter = 2000)
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>'
      message("done")
    }
  )
  return(out)
}
