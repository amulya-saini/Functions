# Required libraries
library(dplyr)
library(Seurat)
library(openxlsx)

compute_gene_exp_cell_prop <- function(
    seurat_object, 
    genes, 
    grouping_column = "Injury", 
    output_file = "compute_gene_exp_cell_prop_summary.xlsx",
    expression_thresholds = list() # This will now accept thresholds for each gene individually
) {
  #' Process and Export Gene Expression Data with gene-specific thresholds
  #'
  #' @param seurat_object A Seurat object containing single-cell RNA-seq data.
  #' @param genes A character vector of gene names to process.
  #' @param grouping_column A metadata column (e.g., "Injury", "Timepoint") to group cells by.
  #' @param output_file A string specifying the name of the Excel file to save the results.
  #' @param expression_thresholds A named list where each gene has its own set of thresholds:
  #'     list(gene_name = c(NoExpression, LowExpression, HighExpression))
  #'     Default thresholds: No Expression = 0, Low = 2, High = 5
  #' @return A named list of pivot tables for each gene.
  
  gene_pivot_tables <- list()
  
  for (gene in genes) {
    # Use the thresholds for each gene, if provided
    if (is.null(expression_thresholds[[gene]])) {
      # If no gene-specific thresholds are provided, use default ones
      thresholds <- c(0, 2, 3)  # Default thresholds: No Expression = 0, Low = 2, High = 3
    } else {
      thresholds <- expression_thresholds[[gene]]
    }
    
    # Extract expression data for the current gene and grouping column
    expression_data <- FetchData(
      seurat_object, 
      vars = c("seurat_clusters", grouping_column, gene)
    )
    
    # Categorize cells based on the current gene's thresholds
    expression_data <- expression_data %>%
      mutate(
        expression_category = case_when(
          get(gene) == expression_thresholds[[gene]][1] ~ "No Expression",
          get(gene) < expression_thresholds[[gene]][2] ~ "Low Expression",
          get(gene) >= expression_thresholds[[gene]][3] ~ "High Expression",
          TRUE ~ "No Expression"  # Default to No Expression for any other case
        )
      )
    
    # Calculate cell counts and proportions for each category
    summary_table <- expression_data %>%
      group_by(seurat_clusters, .data[[grouping_column]], expression_category) %>%
      summarise(cell_count = n(), .groups = "drop") %>%
      group_by(seurat_clusters, .data[[grouping_column]]) %>%
      mutate(proportion = cell_count / sum(cell_count)) %>%
      ungroup()
    
    # Reshape the data for easier interpretation
    pivot_table <- summary_table %>%
      tidyr::pivot_wider(
        names_from = expression_category,
        values_from = c(cell_count, proportion),
        names_sep = "_"
      )
    
    # Add a total cell count column
    pivot_table <- pivot_table %>%
      mutate(
        Total_cell_count = rowSums(
          select(., starts_with("cell_count_")), 
          na.rm = TRUE
        )
      )
    
    # Reorder columns for clarity
    pivot_table <- pivot_table %>%
      select(seurat_clusters, .data[[grouping_column]], Total_cell_count, 
             `cell_count_No Expression`, 
             `cell_count_Low Expression`, 
             `cell_count_High Expression`, 
             everything())
    
    # Store the pivot table in the list
    gene_pivot_tables[[gene]] <- pivot_table
  }
  
  # Write all pivot tables to an Excel file, each as a separate sheet
  wb <- createWorkbook()
  for (gene in names(gene_pivot_tables)) {
    addWorksheet(wb, sheetName = gene)
    writeData(wb, sheet = gene, gene_pivot_tables[[gene]])
  }
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  # Return the list of pivot tables
  return(gene_pivot_tables)
}

# Example Usage
genes_of_interest <- c("gene1", "gene2", "gene3") # Replace with your genes of interest
output_file_path <- "gene_expression_analysis.xlsx" # Path to save the Excel file

# Process the data and save to an Excel file
gene_pivot_tables <- compute_gene_exp_cell_prop(
seurat_object = muscle,
genes = genes_of_interest,
expression_thresholds = list("gene1" = c(0,2,2),
                            "gene2" = c(0,2,2),
                            "gene3" = c(0,2,2)),
grouping_column = "Injury", # Change to "Timepoint" or "Orig.ident" or as needed
output_file = output_file_path
)

# Output: Access the pivot table for a specific gene
View(gene_pivot_tables[["gene1"]])
