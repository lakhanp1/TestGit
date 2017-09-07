library(dplyr)
library(data.table)
library(tibble)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)
library(ggplot2)
require(doParallel)
library(amap)
library(imputeTS)




##################################################################################
## colors for heatmap
colors_seq9 = data.frame(
  YlGnBu = c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58"),
  YlOrRd = c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"),
  RdPu = c("#fff7f3", "#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"),
  OrRd = c("#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000"),
  Reds = c("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d"),
  Oranges = c("#fff5eb", "#fee6ce", "#fdd0a2", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#a63603", "#7f2704"),
  PuBu = c("#fff7fb", "#ece7f2", "#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#045a8d", "#023858"),
  Blues = c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"),
  Purples = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d"),
  GnBu = c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081"),
  stringsAsFactors = F
)


## profile matrix column names
header = c("chr", "start", "end", "geneID", "length", "strand", paste(c("u"), 1:200, sep = ""), paste(c("g"), 1:200, sep = ""), paste(c("d"), 1:100, sep = ""))



## function to impute NA values for the rows with NAs
naImpute = function(index, mat){
  xNew = na.mean(mat[index, ])
  return(xNew)
}

##################################################################################



##################################################################################

## get the profile matrix of class "normalizedMatrix"
getDeeptoolsProfileMatrix = function(file, colNames, signalName){
  
  ## read the profile matrix which is in .gz file
  z1 = paste0("bash -c ", '"zcat ', file, '"', collapse = "")
  
  profileDf = fread(z1, sep = "\t", header = F, skip = 1, na.strings = "nan", col.names = header, data.table = F) %>%
    select(-c(chr, start, end, length, strand)) %>%
    column_to_rownames(var = "geneID")
  
  
  profileMat = data.matrix(profileDf)
  
  ## find rows with NA values and ipmute the NA values: replace NA with mean(row)
  naRows = which(apply(profileMat, 1, function(x) any(is.na(x))))
  profileMat[naRows, ] = do.call(rbind, lapply(naRows, naImpute, mat = profileMat))
  
  ## compare the imputed values with non-imputed values
  # rowN = 13
  # plotNA.imputations(as.numeric(profileDf[naRows[rowN], ]), as.vector(profileMat[naRows[rowN], ]))
  
  ##set attributes for the profileMat to make it of class "normalizedMatrix"
  attr(profileMat, "upstream_index") = 1:200
  attr(profileMat, "target_is_single_point") = FALSE
  attr(profileMat, "target_index") = 201:400
  attr(profileMat, "downstream_index") = 401:500
  attr(profileMat, "extend") = c(2000, 1000)
  attr(profileMat, "signal_name") = signalName
  attr(profileMat, "target_name") = "TSS"
  attr(profileMat, "empty_value") = "NA"
  
  class(profileMat) = c("normalizedMatrix", "matrix")
  
  
  # returns: profile matrix of class "normalizedMatrix"
  return(profileMat)
}

##################################################################################



##################################################################################
## main profile heatmap 
primaryProfileHeatmap = function(profileMat, signalName, numClust, color){
  
  ## hirerchical clustering
  # dend = hclust(dist(profileMat), method = "ward.D")
  # clusterCut = cutree(dend, numClust)
  # clusterData = data.frame(gene = names(clusterCut), cluster = clusterCut, stringsAsFactors = F, row.names = names(clusterCut))
  
  ####### OR ######
  
  ##k-means clustering
  km = kmeans(x = profileMat, centers = numClust, nstart = 100)
  clusterData = data.frame(gene = rownames(profileMat), cluster = NA, stringsAsFactors = F, row.names = rownames(profileMat))
  clusterData[names(km$cluster), "cluster"] = km$cluster
  
  
  clusterData = mutate(clusterData, cluster=paste("Cluster_", cluster, sep = ""))
  row.names(clusterData) = clusterData$gene
  
  # set colors for clusters
  div9_col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#999999", "#ffff33")
  clusterNames = sort(unique(clusterData$cluster))
  clusterColor = structure(div9_col[1:length(clusterNames)], names = clusterNames)
  
  ## annotation legend
  lgd = Legend(at = names(clusterColor), title = "Clusters", 
               type = "lines", legend_gp = gpar(col = clusterColor))
  
  ## generate profile heatmap
  ht = EnrichedHeatmap(
    profileMat, col = color,
    name = signalName, axis_name_rot = 90,
    column_title = signalName,
    split = clusterData["cluster"],
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = clusterColor), show_error = T)),
    top_annotation_height = unit(2, "cm"),
    row_title_rot = 0
  )
  
  ## generate row annotation of clusters
  htAnno = rowAnnotation(df = clusterData["cluster"],
                         col = list(cluster = clusterColor),
                         width = unit(0.5, "cm"),
                         show_legend = F
  )
  
  ##returns: profile heatmap, annotation legend, clusterData
  return(list(
    "heatmap" = ht,
    "rowAnno" = htAnno,
    "legend" = lgd, 
    "cluster" = clusterData,
    "clusterColor" = clusterColor
  ))
}

##################################################################################



##################################################################################
## secondary profile heatmap

secondaryProfileHeatmap = function(profileMat, signalName, clusterData, heatmapColor, clusterColor){
  
  ## annotation legend
  lgd = Legend(at = names(clusterColor), title = "Clusters", 
               type = "lines", legend_gp = gpar(col = clusterColor))
  
  
  ## generate profile heatmap
  ht = EnrichedHeatmap(
    profileMat, col = heatmapColor,
    name = signalName, axis_name_rot = 90,
    column_title = signalName,
    split = clusterData["cluster"],
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = clusterColor), show_error = T)),
    top_annotation_height = unit(2, "cm"),
    row_title_rot = 0
  )
  
  ## generate row annotation of clusters
  htAnno = rowAnnotation(df = clusterData["cluster"],
                         col = list(cluster = clusterColor),
                         width = unit(0.5, "cm"),
                         show_legend = F
  )
  
  ##returns: profile heatmap, annotation legend, clusterData
  return(list(
    "heatmap" = ht,
    "rowAnno" = htAnno,
    "legend" = lgd
  ))
}

##################################################################################



##################################################################################
## read polII expression data
get_secondary_polII_data = function(file, title, clusterData){
  
  polII_df = fread(input = file, header = F, drop = c(1,2,3,5,6), col.names = c("geneID", title),  stringsAsFactors = F, sep = "\t", data.table = F) 
  
  
  clusterData = clusterData %>% left_join(y = polII_df, by = c("gene" = "geneID"))
  
  polII_df = polII_df %>% column_to_rownames(var = "geneID")
  polII_mat = data.matrix(polII_df)
  
  polII_log2_mat = log2(polII_mat + 1)
  
  polII_q = quantile(polII_log2_mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)
  print(paste("Quantiles for polII data", file, sep = ": "), quote = F)
  print(polII_q)
  
  # return: log2(polII expression) matrix, clusterData
  return(list(
    "log2_mat" = polII_log2_mat,
    "clusterDf" = clusterData,
    "quantiles" = polII_q
  ))
}

##################################################################################



##################################################################################
## get the top 10% genes as expressed genes from polII data
get_expressed_polII_genes = function(df, n, field){
  # this function uses interp() from lazyeval package to interpret the dynamic column name for dplyr function top_n()
  topFraction <- round(nrow(df) / n)
  
  call <- quote(df %>%
                  top_n(topFraction, wt = .wt))
  call <- interp(call, .wt = as.name(field))
  newDf <- eval(call)
  return(newDf)
  
}
##################################################################################



##################################################################################
## generate complexHeatmap of single polII data
polII_heatmap = function(log2_matrix, title, color){
  ## column name annotation for allele_ht Heatmap
  colNameAnn <- HeatmapAnnotation(colName = anno_text(x = c(title),
                                                      rot = 90, just = "right",
                                                      offset = unit(1, "npc") - unit(2, "mm")))
  
  ## draw heatmap
  ht = Heatmap(log2_matrix, name = title,
               col = color,
               bottom_annotation = colNameAnn, 
               bottom_annotation_height = unit.c(unit(5, "mm"), max_text_width(title) + unit(2, "mm")),
               show_row_names = FALSE, show_column_names = FALSE,
               width = unit(1, "cm"))
  
  return(ht)
}

##################################################################################



##################################################################################
## add SM_gene information to the expression profile clusterData
add_SM_cluster_info = function(file, clusterDf){
  smData = fread(input = file, header = T, drop = c(2,3,4,5,6,7),  stringsAsFactors = F, sep = "\t", data.table = F)
  smData$is_SM_gene = TRUE
  
  clusterDf = clusterDf %>%
    left_join(y = smData, by = c("gene" = "SM_gene"))
  
  return(clusterDf)
}

##################################################################################




##################################################################################
## add MACS2 peak calling results to expression profile clusterData
add_macs2Peak_info = function(file, clusterDf){
  df = fread(input = file, header = F, select = c(1,2,3,4,7,8,9,13,16),  
             stringsAsFactors = F, sep = "\t", data.table = F,
             col.names = c("peakChr", "peakStart", "peakEnd", "peakId", "peak_enrichment", "peak_pVal", "peak_qVal", "peak_gene", "peak_dist"))
  
  df$isPeak = TRUE
  
  clusterDf = clusterDf %>%
    left_join(y = df, by = c("gene" = "peak_gene"))
  
  return(clusterDf)
  
}

##################################################################################






