# Author: Zheng Li
# Date: 2021-06-05
# Visualizations
#' Plot cell types and each subplot indicates one type
plotCellTypes <- function(xys, labels, ncol)
{
  df <- data.frame(xys, c = labels)
  colnames(df) <- c("x", "y", "c")
  ggplot(df, aes(x, y, color = c)) + 
    geom_point() + 
    facet_wrap(~c, ncol = ncol) +
    gghighlight() +  
    theme_void() + 
    coord_fixed() +
    scale_color_viridis_d("c") 
}


#' concordance matrix between two sets of labels
plotConcordMtx <- function(mtx, methods)
{
  df <- data.frame(
    expand.grid(1:nrow(mtx), 1:ncol(mtx)),
    as.numeric(mtx))
  colnames(df) <- c("x", "y", "counts")
  df$text <- ifelse(df$counts == 0, NA, df$counts)

  ggplot(df, aes(x = y, y = x, fill = counts)) +
    geom_tile(color = "black") +
    geom_text(aes(label = text)) +
    theme_void() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(),
      axis.title.y = element_text(angle = 90),
      axis.text.x = element_text(angle = 90),
      axis.text.y = element_text()
      ) +
    scale_x_discrete(
      name = paste(methods[1], "annotations"),
      limits = colnames(mtx)
      ) +
    scale_y_discrete(
      name = paste(methods[2], "annotations"),
      limits = rownames(mtx)
      ) +
    scale_fill_gradientn(
      values = rescale(c(0, 1, 40, max(df$counts))),
      colors = c("#FFFFFF", "#f7fbff", "#2171b5", "#045a8d"))
}


#' expression of marker genes stratified by cell types
plotMkrs <- function(cnts, markers, labels, legend_title, scale = T)
{
  x <- scater::normalizeCounts(cnts, log = T)
  x <- x[markers, ]
  df <- data.frame(t(x), c = labels)
  df <- df %>% 
    group_by(c) %>%
    summarise_all(mean)
  cnames <- df$c
  df <- df %>% column_to_rownames("c")
  # scale the expression of each gene
  if(scale)
  {
    df <- sapply(df, function(x){
      (x - min(x)) / (max(x) - min(x))
      })
  } else{
    df <- as.matrix(df)
  }
  df <- data.frame(
    expand.grid(1:length(cnames), 1:length(markers)),
    as.numeric(df))
  colnames(df) <- c("x", "y", "Expression")

  ggplot(df, aes(x = y, y = x, fill = Expression)) +
    geom_tile(color = "black") +
    theme_void() +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.text.y = element_text()
      ) +
    scale_x_discrete(limits = markers) +
    scale_y_discrete(limits = cnames) +
    scale_fill_gradient(name = legend_title, low = "#FFFFFF", high = "#e31a1c")
}


#' cell type compositon matrix
plotCellTypeComp <- function(mtx)
{
  df <- data.frame(
    expand.grid(1:nrow(mtx), 1:ncol(mtx)),
    round(as.numeric(mtx), digits = 2)
    )
  colnames(df) <- c("x", "y", "prop")

  ggplot(df, aes(x = y, y = x, fill = prop)) +
    geom_tile(color = "black") +
    geom_text(aes(label = prop)) +
    theme_void() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(),
      axis.text.y = element_text()
      ) +
    scale_x_discrete(limits = colnames(mtx)) +
    scale_y_discrete(limits = rownames(mtx)) +
    scale_fill_gradientn(
      values = scales::rescale(c(0, 0.1, 1)),
      colors = c("#FFFFFF", "#9ecae1", "#2171b5")
      )
}


#' plot cell types or tissue structures
plotClusters <- function(
  xy, 
  labels, 
  title = "", 
  point_size = 1, 
  flip_x = F, 
  flip_y = F, 
  flip_xy = F, 
  ratio = 1
  )
{
  if(flip_xy) xy <- xy[, c(2, 1)]
  if(flip_x) xy[, 1] <- -xy[, 1]
  if(flip_y) xy[, 2] <- -xy[, 2]
  df <- data.frame(
    x = xy[, 1], 
    y = xy[, 2], 
    z = labels
    )
  df$z <- as.factor(df$z)
  
  ggplot(df, aes(x, y, color = z)) + 
    geom_point(size = point_size) + 
    ggtitle(title) + 
    theme_void() + 
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
      ) + 
    coord_fixed(ratio)
}


#' plot ARI boxplots
plotARIBox <- function(BASS_ARI, HMRF_ARI, HMRF_PCA_ARI, ylim = c(0, 1))
{
  cols <- c("#377eb8", "#009E73")
  df <- data.frame(
    ARI = c(HMRF_ARI, HMRF_PCA_ARI),
    Method = c(
      rep("HMRF", length(HMRF_ARI)), 
      rep("HMRF_PCA", length(HMRF_PCA_ARI))
      )
    )

  ggplot(df, aes(x = Method, y = ARI, fill = Method)) + 
    geom_boxplot() + 
    geom_hline(
      yintercept = BASS_ARI, 
      size = 1, 
      linetype = "dashed", 
      color = "red3"
      ) +
    geom_text(
      aes(x = 0.7, y = (BASS_ARI - .04), label = fontface), 
      label = "BASS"
      ) +
    ggtitle("") +
    xlab("") + 
    ylim(ylim) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      legend.position = "none"
      ) +
    scale_fill_manual(values = cols) 
}


#' plot ARI barplot
plotARIBar <- function(ARIs, methods, cols)
{
  df <- data.frame(ARI = ARIs, method = methods)
  ggplot(df, aes(x = method, y = ARI, fill = method)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(
      legend.position = c(0.8, 0.8),
      legend.background = element_blank()) +
    geom_text(aes(label = round(ARI, digits = 2)), vjust = -0.3) +
    scale_fill_manual(values = cols) +
    xlab("") +
    ylim(c(0, 1))
}


#' plot gene expression values
plotGene <- function(xys, cnts, log = T, scale = F, title = "", legend = "", 
  point_size = 1, flip_x = FALSE, flip_y = FALSE, flip_xy = FALSE)
{
  if(flip_xy) xys <- xys[, c(2, 1)]
  if(flip_x) xys[, 1] <- -xys[, 1]
  if(flip_y) xys[, 2] <- -xys[, 2]
  if(log){
    cnts <- log(cnts + 1)
  }
  if(scale){
    cnts <- cnts - min(cnts)
    cnts <- cnts / max(cnts)
  }
  dat_plot <- data.frame(x = xys[, 1], y = xys[, 2], cnts = cnts)
  ggplot(dat_plot, aes(x = x, y = y, color = cnts)) + 
    geom_point(size = point_size) + 
    coord_fixed() + 
    theme_bw() +
    scale_color_gradient(low = "grey60", high = "red") + 
    ggtitle(title) +
    xlab("") +
    ylab("") +
    labs(color = legend) +
    theme(
      plot.title = element_text(hjust = 0.5), 
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())
}


#' plot UMAP
seu_umap <- function(cnts, labels, method, npcs = 50, dims = 1:20)
{
  seu <- CreateSeuratObject(counts = cnts, min.cells = 1)
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu, features = rownames(seu))
  seu <- RunPCA(seu, features = rownames(seu), npcs = npcs)
  seu <- RunUMAP(seu, dims = dims)
  Idents(seu) <- labels
  DimPlot(seu, reduction = "umap", label = T, repel = T) +
    geom_label(aes(x = Inf, y = Inf, hjust = 1, vjust = 1, label = method),
      label.size = NA) +
    theme_classic() +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme(
      legend.position = "None",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.border = element_rect(fill = NA),
      axis.line = element_blank()
      )
}