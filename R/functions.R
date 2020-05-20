#' filter rcorr object and prepare for geomnet plotting
#'
#' @param rcorr.list The object returned by the rcorr() function
#' @param pcut The pvalue cutoff, all correlations with pvalues above this will be excluded
#' @param spearcut The spearman correlation coefficient cutoff.  All correlations with spear values less than this will be removed
#'   Needs improvement as this function cannot only return strong negative correlations...
#'
#' @return Returns an edgelist dataframe, should be basically ready for use with geomnet
#' @export
#'
#' @examples #soon
rcorr_to_geomnet <- function(rcorr.list, pcut=0.05, spearcut=0.6){
  require(reshape2)
  pval <- as.data.frame(rcorr.list$P)
  sim <- as.data.frame(rcorr.list$r)

  sim$to <- rownames(sim)
  pval$to <- rownames(pval)

  pval <- melt(pval, id.vars = 'to')
  sim <- melt(sim, id.vars = 'to')

  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spear')

  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]

  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spear > spearcut,]

  return(sigcor)

}

#' filter ccrepe object and prepare for geomnet plotting
#'
#' @param ccrepe.obj an object returned by the function ccrepe()
#' @param pcut The pvalue cutoff, all correlations with pvalues above this will be excluded
#' @param spearcut The spearman correlation coefficient cutoff.  All correlations with spear values less than this will be removed
#'   Needs improvement as this function cannot only return strong negative correlations...
#' @return Returns an edgelist dataframe, should be basically ready for use with geomnet
#' @export
#'
#' @examples #soon
ccrepe_to_geomnet <- function(ccrepe.obj, pcut=0.05, spearcut=0.6){
  require(reshape2)
  pval <- as.data.frame(ccrepe.obj$p.values)
  sim <- as.data.frame(ccrepe.obj$sim.score)

  sim$to <- rownames(sim)
  pval$to <- rownames(pval)

  pval <- melt(pval, id.vars = 'to', factorsAsStrings=TRUE)
  sim <- melt(sim, id.vars = 'to', factorsAsStrings=TRUE)

  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spear')

  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]

  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spear > spearcut,]
  sigcor$from <- as.character(sigcor$from)
  sigcor$to <- as.character(sigcor$to)
  return(sigcor)

}

#' Extract mothur SILVA formatted taxonomy
#'
#' @param filename A taxonomy file output directly from mothur
#'
#' @return a dataframe with the OTU taxonomy classifications split into all the various taxonomic levels
#' @export
#'
#' @examples #soon
extract_mothur_tax <- function(filename){


  tax <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
  tna <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
  tax[tna] <- NA
  rowcount <- 1

  for (i in strsplit(tax$Taxonomy, split = ';')){
    levelcount <- 1
    for (x in i){
      tax[rowcount,3+levelcount] <- x
      levelcount <- levelcount + 1
      #print(x)
    }
    rowcount <- rowcount+1
  }
  return(tax)
}

#' swap OTU designations with 'genus:OTU'
#'
#' @param mothur_tax A mothur taxonomy file read in by "extract_mothur_tax()"
#'
#' @return Returns a named vector.  The values of this vector are new names for each OTU.  These are made by concatenating the entry in the
#'  'Genus' column with the entry in the 'OTU' column separated by a ':'.  This named vector can then be used to get the new names for
#'  any subset of these OTUs. I use this function when I've got some plot involving a few OTUs and I want to label them with both an OTU number and
#'  some taxonomy information.
#' @export
#'
#' @examples #soon
otu_tax_labels <- function(mothur_tax){
  tax <- mothur_tax
  swip <- tax$OTU
  names(swip) <- paste(tax$Genus, tax$OTU, sep = ':')
  names(swip)[grepl("NA", names(swip))] <- paste(tax$Family[grepl("NA", names(swip))], tax$OTU[grepl("NA", names(swip))], sep = ':')
  swap <- names(swip)
  names(swap) <- swip
  return(swap)
}

#' Remove small unconnected subnetworks from a network
#'
#' @param fortified.edgelist This dataframe should be a network stored as a "fortified.edgelist" (see geomnet package for details?)
#' @param node.dataframe This dataframe should be one returned by gather_nodes() or a bunch of these concatenated together
#' @param min.vert A number specifiing the minimum number of verticies that can be in any subnetwork within the larger network
#'
#' @return Returns a network fit for plotting with geomnet in which small subnetworks with less than the specified number have been removed
#'
#' @export
#'
#' @examples #soon
prune_graph <- function(fortified.edgelist, node.dataframe, min.vert = 1){
  #hopefully this will take a fortified edgelist correlation graph from geomnet, convert it to an igraph object,
  #make sure it's undirected? is this necessary?
  #decompose the graph, throwing away subnetworks that have less than min.vert vertexes
  #recompose the graph and convert it back to the fortified edgelist format
  require(igraph)
  filtered <- graph_from_data_frame(fortified.edgelist, directed = FALSE)
  filtered <- as.undirected(filtered, mode= "collapse")
  filtered <- delete.vertices(filtered, 'NA')
  filtered <- decompose(filtered, mode = "weak", max.comps = NA, min.vertices = min.vert)
  temp <- make_empty_graph(n=0, directed=FALSE)

  for (x in 1:length(filtered)){
    temp <- temp %du% filtered[x]
  }

  filtered <- as_long_data_frame(temp)
  rm(temp)

  filtered$from <- filtered[,length(colnames(filtered))-1]
  filtered$to <- filtered[,length(colnames(filtered))]

  filtered <- filtered[,1:4]
  filtered <- fortify(as.edgedf(filtered), node.dataframe)
  return(filtered)
}

#' Collect type information for nodes in a network
#'
#' @param x a dataframe with columns as measured variables that could be nodes in a correlation network
#' @param typ either a single text value or a text vector of length equal to the number of columns
#'  describing the type(s) of these potential nodes
#'
#' @return returns a two column dataframe. The column 'node' is populated with the column names from the original dataframe.
#'  The 'type' column is user supplied.
#' @export
#'
#' @examples #soon
gather_nodes <- function(x, typ=NA){
  x <- as.data.frame(colnames(x))
  colnames(x)[1] <- 'node'
  x$type <- typ
  return(x)

}

#' this function taken from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
#'
#' @param x An OTU table with taxa as columns and samples as rows
#' @param factors a vector containing the groups you want compared
#' @param sim.method distance metric to use for (dis)similarity values, default is 'bray'
#' @param p.adjust.m method to use for pvalue correction.
#' @param permutations number of permutations to use
#'
#' @return returns a dataframe containing the PERMANOVA rsq, fstat, and pvalue for each pariwise comparison
#' @export
#'
#' @examples #soon
pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m = 'none', permutations=999){
  #this function taken from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R

  require(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()

  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = permutations);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


#' Generate NMDS points, group centroids and standard error ellipses
#'
#' @param metadata A dataframe containing the metadata relating to the samples in your OTU table.  Seems to break if this dataframe only has 1 column.
#' @param OTU_table an OTU table with samples as rows and taxa as columns
#' @param grouping_set a column name from your metadata.  This column should contain information
#'   about how your samples should be grouped for centroid and ellipse calculation
#' @param distance_method a single text value that describes the desired distance metric to use, must be
#'  accepted by the function vegdist()
#' @param rand_seed numeric value to use as a random seed (for reproducability?)
#' @param MDS_trymax number of iterations to use for the metaMDS call, default is 1000
#' @param autotransform passed to vegan's metaMDS, should community data be transformed?
#'     Should be TRUE/FALSE.  Default = FALSE
#' @param wascores passed to vegan's metaMDS, should species scores be calculated?
#'     Should be TRUE/FALSE.  Default = TRUE
#' @param expand passed to vegan's metaMDS, should species scores be expanded?
#'
#'
#' @return Returns a list containing 3 items. [[1]] = input metadata with NMDS coordinates, group centroid coordinates for each sample.
#'   [[2]] = a dataframe containing information about the standard error confidence ellipses for each group.  This datafame
#'   contains 100 points per group and is suitable for using with geom_path().
#'   [[3]] = the ordination object created by the vegan function metaMDS().
#'
#' @export
#'
#' @examples #soon
NMDS_ellipse <- function(metadata, OTU_table, grouping_set,
                         distance_method = 'bray',
                         rand_seed = 77777,
                         MDS_trymax = 1000,
                         autotransform = FALSE,
                         wascores = TRUE,
                         expand = FALSE){
  require(vegan)
  require(tidyr)

  if (grouping_set %in% colnames(metadata)){
    if (all(rownames(metadata) == rownames(OTU_table))){

      set.seed(rand_seed)
      generic_MDS <- metaMDS(OTU_table, k = 2,
                             trymax = MDS_trymax,
                             autotransform = autotransform,
                             distance = distance_method,
                             wascores = wascores,
                             expand = expand)

      stress <- generic_MDS$stress
      nmds_points <- as.data.frame(generic_MDS$points)
      metadata <- metadata[match(rownames(generic_MDS$points), rownames(metadata)),]
      metadata <- as.data.frame(metadata) # weird things were happening when a grouped tibble was used as metadata...
      metanmds <- cbind(metadata, nmds_points)
      # browser()
      nmds.mean <- aggregate(metanmds[,grep("MDS", colnames(metanmds))], list(group=metanmds[[grouping_set]]), mean)
      metanmds[[grouping_set]] <- factor(metanmds[[grouping_set]]) # this 'set' needs to be passed in from function

      #check to make sure at least 3 obs for each grouping_set

      numobs <- metanmds %>% group_by(!!grouping_set) %>% summarise(n=n())
      if (min(numobs$n) >= 3){
        ord <- ordiellipse(generic_MDS, metanmds[[grouping_set]], label = TRUE, conf = .95, kind = 'se', draw = 'none')

        df_ell <- data.frame()
        for (d in levels(metanmds[[grouping_set]])){
          df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds[[grouping_set]] == d,],
                                                           vegan:::veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
        }



        # this loop assigns the group centroid X coordinates to each sample, there is probably a better way...

        metanmds$centroidX <- NA
        metanmds$centroidY <- NA


        for (level in levels(metanmds[[grouping_set]])){
          metanmds[metanmds[[grouping_set]] == level,]$centroidX <- nmds.mean$MDS1[nmds.mean$group == level]
          metanmds[metanmds[[grouping_set]] == level,]$centroidY <- nmds.mean$MDS2[nmds.mean$group == level]


        }
        print(paste('Ordination stress:', stress, sep = ' '))
        return(list(metanmds, df_ell, generic_MDS))

      } else {
        warning('One of your groups in "grouping_set" has less than 3 observations, cannot generate elipses')
        df_ell <- data.frame()
        return(list(metanmds, df_ell, generic_MDS))}

    } else {
      stop('The rownames for your OTU table and metadata do not match.')
    }

  }else {
    stop('The grouping set column you have provided in not in your metadata.')
  }



}



#' Quickly plot Deseq2 results
#'
#' Should really split this into 2 functions, one that returns the results df with lfc shrinkage etc
#'  and another that plots...
#'
#' @param DeSeq.object Deseq object that the desired pairwise contrast exists in
#' @param phyloseq.object phyloseq object containing the taxonomy table associated with the OTUs in the Deseq object
#' @param pvalue pvalue to filter results table with
#' @param taxlabel the taxonomy label to use for the plot labels (needs to be present in the tax table from the phyloseq object)
#' @param colors optional. A vector of length 2, first value will be the color associated with the 2nd value from the contrast vector, 2nd value will be the color associated with the 3rd value from the contrast vector
#' @param name the name of the results contrast you want to display
#' @param shrink_type one of the methods used by the DEseq function lfcshrink
#' @param scientific do you want pvalues in scientific format?
#' @param cookscut should we apply the cookscutoff parameter in DEseq2?
#' @param alpha passed to DEseq2 results function
#' @param LFC_cut filter results based on absolute value lfc? default=0
#'
#' @return returns a list containing: [[1]] a ggplot object, [[2]] a dataframe containing the significantly differentially abundant features
#' @export
#'
#' @examples #soon
Deseq.quickplot <- function(DeSeq.object,phyloseq.object,
                            pvalue = 0.05, name,
                            taxlabel = 'Genus', shrink_type = 'normal' ,
                            colors=NULL, scientific=FALSE,
                            cookscut=FALSE, alpha=0.05,
                            LFC_cut=0){
  require(ggplot2)
  require(phyloseq)
  require(DESeq2)
  tmp1 <- sub('.*_(.*)_vs_(.*)', '\\1', name)
  tmp2 <- sub('.*_(.*)_vs_(.*)', '\\2', name)
  # contrast.vector <- c(tmp1, tmp2, tmp3)
  res <- results(object = DeSeq.object, name = name, cooksCutoff = cookscut, alpha = alpha)
  res <- lfcShrink(DeSeq.object,res=res, coef = name, type = shrink_type)
  sigtab = res[which(res$padj < pvalue), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq.object)[rownames(sigtab), ], "matrix"))
  sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = scientific)
  sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, tmp1, tmp2)
  sigtab$OTU <- rownames(sigtab)
  sigtab <- sigtab[abs(sigtab$log2FoldChange) > LFC_cut,]

  if (is.null(colors)){

    deseq <- ggplot(sigtab, aes(x=reorder(rownames(sigtab), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes_string(x='OTU', y=0, label = taxlabel), size=3)+
      theme(axis.text.x=element_text(color = 'black', size = 10),
            axis.text.y=element_text(color = 'black', size=10),
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ coord_flip() + xlab(element_blank())

  }


  if (!is.null(colors)){
    fofill <- c(colors[1], colors[2])
    names(fofill) <- c(contrast.vector[2], contrast.vector[3])
    deseq <- ggplot(sigtab, aes(x=reorder(rownames(sigtab), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes_string(x='OTU', y=0, label = taxlabel), size=3)+
      scale_fill_manual(values = fofill) +
      theme(axis.text.x=element_text(color = 'black', size = 10),
            axis.text.y=element_text(color = 'black', size=10),
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ coord_flip() + xlab(element_blank())

  }
  return(list(deseq, sigtab))

}


#' CTs to log2 fold change relative to controls
#'
#' @param data this is the input dataframe, should have genes as columns and samples as rows, replicates should be collpsed and only one housekeeping gene
#'  is supported.  I don't think NA's or missing data is supported...not really sure.
#' @param exp_genes_indicies A numeric vector containing the column indicies for the genes of experimental interest
#' @param housekeeping_gene_index A single numeric value defining the column index for the housekeeping gene
#' @param control_row_indicies A numeric vector containing the row indicies for the samples that are in the control group.  Final data will be
#'  expressed relative to the mean of this group for each gene.
#' @param metadata_column_indicies Numeric column indicies defining which columns are metadata and should be returned along with the final results
#'
#' @return returns a dataframe of log2(2^DDct) values centered on the mean of the control group
#' @export
#'
#' @examples #soon
qpcr_DDCT <- function(data, exp_genes_indicies, housekeeping_gene_index, control_row_indicies, metadata_column_indicies){
  delta.CT <- data[,exp_genes_indicies] - data[,housekeeping_gene_index]
  #delta.CT <- delta.CT[,-housekeeping_gene_index]
  gene_av_controls <- colSums(delta.CT[control_row_indicies,])/length(control_row_indicies)  # finds the average delta CT for each gene of the control group only
  delta_delta_CT <- scale(delta.CT, center=gene_av_controls, scale = FALSE) # this is the calculation of the second delta CT, subracts the average control CT for each gene from each respective gene
  delta_delta_CT <- -(delta_delta_CT)
  result <- cbind(delta_delta_CT, data[,metadata_column_indicies])
  return(result)

}

#' CTs to log2 fold change relative to controls with LOD estimates
#'
#' @param data is is the input dataframe, should have genes as columns and samples as rows, replicates should be collpsed and only one housekeeping gene
#'  is supported.  I don't think NA's or missing data is supported...not really sure.
#' @param exp_genes_indicies A numeric vector containing the column indicies for the genes of experimental interest
#' @param housekeeping_gene_index A single numeric value defining the column index for the housekeeping gene
#' @param control_row_indicies A numeric vector containing the row indicies for the samples that are in the control group.  Final data will be
#'  expressed relative to the mean of this group for each gene.
#' @param metadata_column_indicies Numeric column indicies defining which columns are metadata and should be returned along with the final results
#' @param LOD Your desired Limit of detection, all cts above this value will be replaced with NA and then this value
#' @param LOD_conf The confidence level for the LOD estimate
#'
#' @return returns list of length 2, 1) a dataframe of log2(2^DDct) values centered on the mean of the control group
#' and 2) a dataframe of LOD estimates for each gene
#' @export
#'
#' @examples #soon
qpcr_DDCT_LOD <- function(data,
                          exp_genes_indicies,
                          housekeeping_gene_index,
                          control_row_indicies,
                          metadata_column_indicies,
                          LOD, LOD_conf=.999){

  delta.CT <- data[,exp_genes_indicies] - data[,housekeeping_gene_index]

  sampleLODs <- LOD - data[,housekeeping_gene_index]
  #delta.CT <- delta.CT[,-housekeeping_gene_index]
  gene_av_controls <- colSums(delta.CT[control_row_indicies,])/length(control_row_indicies)  # finds the average delta CT for each gene of the control group only

  geneLODs <- LOD - gene_av_controls
  resu <- matrix(ncol=length(geneLODs), nrow = 0)
  colnames(resu) <- names(geneLODs)
  for (i in 1:length(sampleLODs)){
    resu <- rbind(sampleLODs[i] - geneLODs, resu)
  }
  LOD_est <- apply(resu, 2, t.test, conf.level = LOD_conf)
  lowers <- c()
  uppers <- c()

  for (i in 1:length(LOD_est)){
    lowers[i] <- LOD_est[[i]]$conf.int[1]
    uppers[i] <- LOD_est[[i]]$conf.int[2]
  }

  LOD_df <- data.frame(gene=names(LOD_est),
                       l_conf=lowers,
                       u_conf=uppers)


  delta_delta_CT <- scale(delta.CT, center=gene_av_controls, scale = FALSE) # this is the calculation of the second delta CT, subracts the average control CT for each gene from each respective gene
  log2expDDCT <- log(2^(-(delta_delta_CT)), base = 2)  # changes to log *technically this does nothing unless you change the log base
  result <- cbind(log2expDDCT, data[,metadata_column_indicies])

  return(list(result, LOD_df))
}

#' Return a darker hue of a color
#'
#'I did not write this function, taken from:
#' https://gist.github.com/Jfortin1/72ef064469d1703c6b30
#'
#' @param color input color to change
#' @param factor how much darkening to do
#'
#' @return
#' @export
#'
#' @examples #soon
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


#' Return a lighter hue of a color
#'
#'I did not write this function, taken from:
#' https://gist.github.com/Jfortin1/72ef064469d1703c6b30
#'
#' @param color color to lighten
#' @param factor how much lightening to do
#'
#' @return
#' @export
#'
#' @examples #soon
lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}

#' DESeq2 differential abundance for all treatments compared to 'Control'
#'
#' @param phyloseq phyloseq object
#' @param day day of interest (used to subset phyloseq object)
#' @param tissue tissue of interest (used to subset phyloseq object)
#' @param scientific change pvals to sci notation?
#' @param shrink_type what type of lfc shrinkage? passed to DESeq2's lfc_shrink.
#' should be one of 'normal', 'ashr', 'apeglm'.
#' @param alpha passed to DESeq2 results() function.
#' @param cooks_cut apply cooks cutoff?
#' @param pAdjustMethod what type of p.adjustment to apply?
#'
#' @return returns a dataframe containing all the significantly differentially abundant taxa for
#' each treatment in your experiment compared to the 'Controls'
#'
#' @export
#'
#' @examples #soon
DESeq_difabund <- function(phyloseq, day, tissue, scientific = TRUE, shrink_type='normal',
                           alpha=0.1, cooks_cut=FALSE, pAdjustMethod='BH'){

  # FS12b.glom <- tax_glom(FS12b, taxrank = 'Genus')
  FS12b.glom <- prune_samples(x = phyloseq, samples = phyloseq@sam_data$day == day & phyloseq@sam_data$tissue == tissue)
  FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
  FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~treatment)
  FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')

  finres <- list()
  resind <- 1
  for (i in 2:length(resultsNames(FS12.de))){
    print(resultsNames(FS12.de)[i])
    treat <- sub('treatment_(.*)_vs_Control','\\1',resultsNames(FS12.de)[i])
    comp <- sub('treatment_', '', resultsNames(FS12.de)[i])
    res <- results(object = FS12.de, name = resultsNames(FS12.de)[i], alpha=alpha, cooksCutoff = cooks_cut, pAdjustMethod = pAdjustMethod)
    res <- lfcShrink(FS12.de, coef = resultsNames(FS12.de)[i], type = shrink_type)
    sigtab = res[which(res$padj < alpha), ]

    if (nrow(sigtab) != 0){
      # browser()
      sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(FS12b.glom)[rownames(sigtab), ], "matrix"))
      sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = scientific)
      sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, treat, paste('down',treat, sep = '_'))
      sigtab$OTU <- rownames(sigtab)
      sigtab$tissue <- tissue
      sigtab$day <- day
      sigtab$comp <- comp
      finres[[resind]] <- sigtab

      resind <- resind + 1
    }



  }

  finres <- bind_rows(finres)
  return(finres)

}
