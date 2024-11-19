
params_to_df <- function(params){
  ls <- params
  class(ls) <- "list"
  for(i in names(ls)){
    if(length(ls[[i]]) > 1){
      ls[[i]] <- paste(ls[[i]], collapse = "; ")
    } else if(is.character(ls[[i]])){
      if(file.exists(ls[[i]])){
        ls[[i]] <- basename(ls[[i]])
      }
    }
  }
  #return(ls)
  df <- data.frame("parameter" = names(ls),
                   "argument" = unlist(ls, use.names = F)) %>%
    dplyr::filter(argument != "")
  return(knitr::kable(df))
}


DotPlot_Plus <- function(seurat, 
                         features,
                         #slot = "data",
                         assay = NULL, 
                         col.min = -2.5,
                         col.max = 2.5, 
                         dot.min = 0,
                         dot.scale = 6,
                         idents = NULL,
                         group.by = NULL,
                         scale = T,
                         scale.by = "radius",
                         scale.min = NA,
                         scale.max = NA){
  ### inpsired from https://davemcg.github.io/post/lets-plot-scrna-dotplots/
  require(patchwork)
  require(ggtree)
  require(aplot)
  require(viridis)
  require(tidyverse)
  
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  df <- DotPlot(seurat, features = features, group.by = group.by, col.min = col.min, col.max = col.max, dot.min = dot.min, assay = assay, scale = scale)$data %>%
    dplyr::rename(`Percent Expressed` = "pct.exp",
                  `Average Expression` = "avg.exp.scaled")  %>%
    mutate(`Percent Expressed` = ifelse(`Percent Expressed` < 0.05, 0, `Percent Expressed`)) %>%
    mutate(`Average Expression` = ifelse(is.na(`Average Expression`), min(`Average Expression`, na.rm = T), `Average Expression`))
  
  if("feature.groups" %in% colnames(df)){
    df <- df %>%
      mutate(`feature.groups` = as.character(`feature.groups`)) %>%
      mutate(`features.plot` = `feature.groups`)
  }
  
  hclust_prep <- df %>%
    dplyr::select(features.plot, id, `Average Expression`) %>%
    spread(id, `Average Expression`) %>%
    column_to_rownames("features.plot") %>%
    as.matrix()
  
  hclust_rows <- hclust(dist(hclust_prep))
  hclust_cols <- hclust(dist(t(hclust_prep)))
  
  ggtree_rows <- ggtree::ggtree(as.dendrogram(hclust_rows), branch.length = "none")
  ggtree_cols <- ggtree::ggtree(as.dendrogram(hclust_cols)) +  layout_dendrogram()
  
  
  df <- df %>%
    mutate(features.plot = factor(features.plot, levels = hclust_rows$labels[hclust_rows$order])) %>%
    mutate(id = factor(id, levels = hclust_cols$labels[hclust_cols$order]))
  
  p <- ggplot(df, aes(x= id, y = features.plot, color = `Average Expression`, size = `Percent Expressed`)) + 
    geom_point() + 
    cowplot::theme_cowplot() + 
    #theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("Identities") +
    ylab('') +
    #theme(axis.ticks = element_blank()) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    scale_color_gradientn(colours = viridis::viridis(20), 
                          limits = c(floor(min(df$`Average Expression`)),
                                     ceiling(max(df$`Average Expression`))), 
                          oob = scales::squish, name = 'Average Expression') +
    scale_y_discrete(position = "right")
  
  ggtree_rows <- ggtree_rows + ylim2(p)
  ggtree_cols <- ggtree_cols + xlim2(p)
  
  
  p_out <- plot_spacer() +  plot_spacer() + ggtree_cols  +
    plot_spacer() + plot_spacer() + plot_spacer() +
    ggtree_rows +  plot_spacer() + p +
    plot_layout(ncol = 3, widths = c(0.3, -0.1, 4), heights = c(0.3, -0.1,  4))
  return(p_out)
}

subchunkify <- function(g, fig_height=10, fig_width=10) {
  ## @nova https://stackoverflow.com/questions/15365829/dynamic-height-and-width-for-knitr-plots
  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')
  
  sub_chunk <- paste0("`","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=",
                      fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
                      "\n(", 
                      g_deparsed
                      , ")()",
                      "\n`","``")
  
  cat(trimws(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE)))
}


clustree_arrow <- function(clustree_plot, selected_res){
  require(clustree)
  df <- clustree_plot$data[clustree_plot$data[,4] == selected_res,]
  y_sel <- unique(df$y)[1]
  x_sel1 <- max(df$x) + 0.5
  x_sel2 <- max(df$x) + 2
  x_label <- x_sel2 + 0.5
  p <- clustree_plot +
    geom_segment(aes(x = x_sel2 , y = y_sel, xend = x_sel1 , yend = y_sel),
                 arrow = arrow(length = unit(0.3, "cm"))) +
    geom_text(x=x_label, y=y_sel, label=selected_res) 
  return(p)
}

colorAssign <- function(valueVector, scale_limits = NULL, colors = c("blue", "white", "red"), length.vector = 500){
  require(RColorBrewer)
  colorLS <- colorRampPalette(colors = colors)(length.vector)  
  if(is.null(scale_limits)){
    breaks <- seq(-max(abs(valueVector)), max(abs(valueVector)), length.out = length.vector)
  } else {
    breaks <- seq(scale_limits[1], scale_limits[2], length.out = length.vector)
  }
  
  
  minBreak <- which(abs(breaks - min(valueVector)) == min(abs(breaks - min(valueVector))))
  maxBreak <- which(abs(breaks - max(valueVector)) == min(abs(breaks - max(valueVector))))
  return(colorLS[minBreak:maxBreak])
}


FeaturePlot_bidirectional <- function(object, feature, pt.size = NULL, order = FALSE,
                                      reduction = NULL, split.by = NULL, slot = "data",
                                      cols = c("darkblue", "blue", "lightgray", "red", "darkred")) {
  require(ggplot2)
  require(Seurat)
  active_assay <- object@active.assay
  if(!feature %in% row.names(object) ){
    stop(paste("ERROR: Feature `", feature,"` is not in the ", active_assay, "data.", sep = "" ))
    
  }
  colLS <- colorAssign(GetAssayData(object = object, slot = slot)[feature,],  
                       length.vector = 100,
                       colors = cols)
  p <- FeaturePlot(object, features = feature, reduction = reduction, slot = slot, order = order, pt.size = pt.size) 
  p$data <- p$data %>%
    rownames_to_column("cell") %>%
    inner_join(., seurat@meta.data %>%
                 rownames_to_column("cell"), by = "cell") %>%
    column_to_rownames("cell")
  
  p <- p + 
    scale_color_gradientn(colors = colLS) 
  if(!is.null(split.by)){
    if(split.by %in% names(object@meta.data)) {
      p <- p + 
        facet_wrap(as.formula(paste("~", split.by)))
    } else if(!split.by %in% names(object@meta.data)  ){
      stop(paste("ERROR: split.by variable `", split.by,"` is not in the metadata.", sep = "" ))
    }
  }
  return(p)
}

compute_ident_frequencies <- function(seurat, 
                                      ident = "seurat_clusters",
                                      outcomes = NULL){
  require(checkmate)
  require(tidyverse)
  coll <- checkmate::makeAssertCollection()
  if (methods::is(seurat) != "Seurat") {
    coll$push("Seurat object not recognized")
  }
  
  checkmate::assertCharacter(ident, null.ok = F ,.var.name = "ident", add = coll)
  checkmate::assertCharacter(outcomes, null.ok = T, any.missing = F ,.var.name = "outcomes", add = coll)
  
  checkmate::reportAssertions(coll)

  
  df <- seurat@meta.data %>%
    dplyr::select(dplyr::all_of(c(ident, "SampleID"))) %>%
    dplyr::group_by_at(vars(dplyr::all_of(ident),"SampleID" )) %>%
    dplyr::summarise(absolute_counts = dplyr::n()) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(total_n = sum(absolute_counts)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(relative_frequency = absolute_counts/ total_n * 100)
  
  

  if(!is.null(outcomes)){
    var_setdiff <- setdiff(outcomes, colnames(seurat@meta.data))
    if(length(var_setdiff)> 0){
      coll$push(paste0("Error: Sample metadata does not contain variable(s) ", paste(sQuote(var_setdiff), collapse = ";")))
    } else {
      sample_meta <- seurat@meta.data %>%
        dplyr::select(dplyr::all_of(c(outcomes, "SampleID"))) %>%
        unique()
      df <- df %>%
        left_join(sample_meta, by = "SampleID")
    }
    
  }
  checkmate::reportAssertions(coll)
  
  return(df)
  
  
}

cluster_frequency_plot <- function(seurat, 
                                   ident = "seurat_clusters",
                                   split.by = NULL,
                                   ncols = NULL){

  coll <- checkmate::makeAssertCollection()
  if (methods::is(seurat) != "Seurat") {
    coll$push("Seurat object not recognized")
  }
  
  checkmate::assertCharacter(split.by, null.ok = T, any.missing = F ,max.len = 3, .var.name = "split.by", add = coll)
  checkmate::assertNumeric(ncols, null.ok = T, any.missing = F ,len = 1, .var.name = "ncols", add = coll)
  checkmate::reportAssertions(coll)
  
  
  df <- compute_ident_frequencies(seurat,
                                  ident = ident,
                                  outcomes = split.by)
  
  col_extract <- ggplot_build(DimPlot(seurat, group.by = ident))$data[[1]] %>%
                  dplyr::select(colour, group) %>%
                  unique() 
  
  
  if( ident == "seurat_clusters"){
    col_extract$label <- as.character(col_extract$group - 1)
  } else if(class(seurat@meta.data[,ident]) == "character"){
    label_df <- data.frame(label = levels(factor(seurat@meta.data[,ident])),
                           group = seq(1, length(unique(seurat@meta.data[,ident]))))
    col_extract <- col_extract %>%
                      left_join(label_df, by = "group")
  } else if(class(seurat@meta.data[,ident]) == "factor"){
    label_df <- data.frame(label = levels(seurat@meta.data[,ident]),
                           group = seq(1, length(unique(seurat@meta.data[,ident]))))
    col_extract <- col_extract %>%
                left_join(label_df, by = "group")
  }
  
  #return(col_extract)
  colorLS <- col_extract$colour
  names(colorLS) <- col_extract$label
  

  p <-   ggplot(df, aes_string(x = "SampleID", y = "relative_frequency", fill = ident)) +
    theme_classic()  +
    geom_bar(position = "stack", stat = "identity") +
    ylab("Relative frequency (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    scale_fill_manual(values = colorLS)
  
  
  
  if(!is.null(split.by)){
    if(length(split.by) == 1){
      p <- p + facet_wrap(stats::as.formula(paste("~ `", split.by, "`", sep = "")), ncol = ncols, scales = "free" )
    } else {
      split_2p <- paste(unlist(lapply(split.by[c(2:length(split.by))], function(x){
        return(paste("`", x, "`", sep = ""))
      })), collapse = " + ")
      p <- p + facet_grid(stats::as.formula(paste("`", split.by[1], "`", " ~ ", split_2p, sep = "")), scales = "free")
    }
  }
  return(p)
}


