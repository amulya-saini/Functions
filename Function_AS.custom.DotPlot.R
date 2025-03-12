library(cowplot)
library(ggplot2)
library(Seurat)
library(tibble)
library(tidyr)

#' Custom Dot Plot for Single-Cell RNA-Seq Data
#'
#' A custom function to generate a dot plot for visualizing gene expression across different conditions,
#' with options to scale the expression values and customize the plot appearance.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param features A feature (gene) name to be visualized in the dot plot. (accepts only one gene)
#' @param assay The name of the assay to use for the plot (default is `NULL`, in which case the default assay is used).
#' @param split.column The metadata column used to split the data (default is `"orig.ident"`).
#' @param cols A vector of two colors to represent the expression scale (default is `c("lightgrey", "blue")`).
#' @param col.min The minimum value for color scaling (default is `-2.5`).
#' @param col.max The maximum value for color scaling (default is `2.5`).
#' @param dot.min The minimum percentage of cells expressed for dot size scaling (default is `0`).
#' @param dot.scale The scaling factor for dot size (default is `6`).
#' @param scale.by The type of scaling to use (`'size'` or `'radius'` for the dot size scaling, default is `'radius'`).
#' @param scale Logical, whether to scale the expression data (default is `TRUE`).
#' @param scale.min The minimum value for scaling the expression data (default is `NA`).
#' @param scale.max The maximum value for scaling the expression data (default is `NA`).
#'
#' @return A `ggplot` object representing the custom dot plot.
#'
#' @details This function generates a custom dot plot based on a Seurat object, where:
#' - Each dot represents the average expression of a given gene (feature) for cells in a particular condition (split.column),
#' - The dot size corresponds to the percentage of cells expressing the feature (above a given threshold),
#' - The color of the dots represents the average expression level of the feature.
#'
#' The function scales the expression data, reshapes it to a long format for ggplot, and adjusts the color scale based on user-defined limits.
#'
#' @examples
#' # Create a custom dot plot for a Seurat object with specific features and conditions
#' CustomDotPlot(object = seurat_object, 
#'               features = "geneA", 
#'               split.column = "condition", 
#'               cols = c("lightgrey", "blue"))
#'

AS.custom.DotPlot <- function(
    object,
    features,
    assay = NULL,
    split.column = "orig.ident",
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    scale.by = 'radius',
    scale = TRUE,
    scale.min = NA,
    scale.max = NA
) {
  # Set the assay
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  
  # Ensure correct scale function
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  
  # Extract relevant data
  cells <- colnames(object[[assay]])
  data.features <- FetchData(object = object, vars = features, cells = cells)
  
  # Extract split information
  split.values <- object@meta.data[, split.column, drop = TRUE]
  data.features$id <- Idents(object = object)[cells]
  
  if (!is.factor(data.features$id)) {
    data.features$id <- factor(data.features$id)
  }
  
  # Create a combined identifier for each cell
  data.features$split <- split.values
  data.features$id <- paste(data.features$id, data.features$split, sep = "--")
  
  # Compute avg.exp and pct.exp
  data.plot <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident, 1:(ncol(data.features) - 2), drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    pct.exp <- apply(data.use, 2, function(x) PercentAbove(x, threshold = 0))
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(data.plot) <- unique(data.features$id)
  
  # Convert list to dataframe
  data.plot <- do.call(rbind, lapply(names(data.plot), function(x) {
    df <- as.data.frame(data.plot[[x]])
    df$features.plot <- rownames(df)
    df$id <- x
    return(df)
  }))
  
  # Extract split categories for restructuring
  data.plot$split <- sub(".*--", "", data.plot$id)
  data.plot$id <- sub("-.*", "", data.plot$id)
  
  # Reshape to wide format
  data.wide <- reshape(data.plot, idvar = c("features.plot", "id"), timevar = "split", direction = "wide")
  data.wide[is.na(data.wide)] <- 0
  
  # Scale the data
  if (scale) {
    for (split in unique(split.values)) {
      # Ensure column names are properly constructed
      avg_exp_col <- paste0("avg.exp.", split)
      
      # Check if the column exists to prevent errors
      if (avg_exp_col %in% colnames(data.wide)) {
        data.wide[[avg_exp_col]] <- scale(log1p(data.wide[[avg_exp_col]]))
        data.wide[[avg_exp_col]] <- MinMax(data.wide[[avg_exp_col]], min = col.min, max = col.max)
      } else {
        warning(paste("Column", avg_exp_col, "does not exist in data.wide"))
      }
    }
  }
  data.wide[is.na(data.wide)] <- 0
  
  # Reshape to long format for ggplot
  data.long <- pivot_longer(data.wide, 
                            cols = starts_with("avg.exp."), 
                            names_to = "condition", 
                            values_to = "avg.exp")
  
  # Fix pct.exp - Make sure to correctly align pct.exp for each condition
  pct.long <- pivot_longer(data.wide, 
                           cols = starts_with("pct.exp."), 
                           names_to = "pct_condition", 
                           values_to = "pct.exp")
  
  # Clean column names
  data.long$condition <- gsub("avg.exp.", "", data.long$condition)
  pct.long$pct_condition <- gsub("pct.exp.", "", pct.long$pct_condition)
  
  # Ensure pct.exp values are correctly scaled
  pct.long$pct.exp <- pmax(pct.long$pct.exp, dot.min) * 100 
  
  # Merge avg.exp and pct.exp data
  data.long <- merge(data.long, pct.long, by.x = c("features.plot", "id", "condition"), by.y = c("features.plot", "id", "pct_condition"), all = TRUE)
  
  # Set condition (split.column values) as a factor and order based on levels from the metadata
  split.levels <- factor(split.values, levels = levels(split.values))
  data.long$condition <- factor(data.long$condition, levels = levels(split.levels))
  
  # Generate the plot
  plot <- ggplot(data = data.long, mapping = aes(x = condition, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp)) +
    scale_size_continuous(labels = scales::percent) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = split.column,
      y = 'Identity'
    ) +
    theme_cowplot() +
    scale_color_gradient(low = cols[1], high = cols[2]) +
    guides(color = guide_colorbar(title = 'Average Expression'))
  
  return(plot)
}


# Example Usage
# AS.custom.DotPlot(seurat_obj, features = "vcanb", split.column = "orig.ident",
#               scale.max = 100, scale.min = 0,
#               dot.scale = 8)+
#   RotatedAxis() +
#  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
# theme(plot.title = element_text(hjust = 0.5))
