NMDS_ellipse2 <- function(metadata, OTU_table, grouping_set,
                         distance_method = 'bray',
                         rand_seed = 77777,
                         MDS_trymax = 1000,
                         autotransform = FALSE,
                         wascores = TRUE,
                         expand = FALSE,
                         dist_matrix){
  require(vegan)
  if (grouping_set %in% colnames(metadata)){
    if (exists(dist_matrix)){
      if (class(dist_matrix) == 'dist'){
        if (all(attributes(dist_matrix)$Labels == rownames(metadata))){
          set.seed(rand_seed)
          generic_MDS <- metaMDS(dist_matrix, k = 2,
                                 trymax = MDS_trymax)

        } else {
          stop('Metadata rownames do not match rownames of provided distance matrix')
        }
      } else {
        stop('Your input distance matrix is not of class "dist"')
      }
    }
    if (exists(dist_matrix) & exists(OTU_table)){
      stop('You have provided both a distance matrix and an OTU table, please provide only one or the other')
    }
    if (exists(OTU_table)){

    }
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
      metanmds <- cbind(metadata, nmds_points)
      nmds.mean <- aggregate(metanmds[,grep("MDS", colnames(metanmds))], list(group=metanmds[[grouping_set]]), mean)
      metanmds[[grouping_set]] <- factor(metanmds[[grouping_set]]) # this 'set' needs to be passed in from function

      #check to make sure at least 3 obs for each grouping_set

      numobs <- metanmds %>% group_by_(grouping_set) %>% summarise(n=n())
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

  } else {
    stop('The grouping set column you have provided in not in your metadata.')
  }



}
