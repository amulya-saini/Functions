# Required libraries
library(dplyr)
library(Seurat)
library(openxlsx)

# Function to process gene expression and save to an Excel file
compute_gene.exp_cell.prop <- function(
    seurat_object, 
    genes, 
    grouping_column = "Injury", 
    output_file = "gene_expression_summary.xlsx"
    ) {
  #' Process and Export Gene Expression Data
  #'
  #' This function processes gene expression data from a Seurat object, categorizes cells 
  #' by expression levels for each gene of interest, and writes the results to an Excel file.
  #' Each gene's results are saved in a separate sheet within the Excel file.
  #'
  #' @param seurat_object A Seurat object containing single-cell RNA-seq data.
  #' @param genes A character vector of gene names to process.
  #' @param grouping_column A metadata column (e.g., "Injury", "Timepoint") to group cells by.
  #' @param output_file A string specifying the name of the Excel file to save the results.
  #'
  #' @return A named list of pivot tables for each gene.
  
  # Initialize a list to store pivot tables for each gene
  gene_pivot_tables <- list()
  
  # Loop through each gene
  for (gene in genes) {
    # Extract expression data for the current gene and grouping column
    expression_data <- FetchData(
      seurat_object, 
      vars = c("seurat_clusters", grouping_column, gene)
    )
    
    # Categorize cells based on gene expression levels
    expression_data <- expression_data %>%
      mutate(
        expression_category = case_when(
          get(gene) == 0 ~ "No Expression",
          get(gene) > 2 ~ "High Expression",
          TRUE ~ "Low Expression"
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
gene_pivot_tables <- compute_gene.exp_cell.prop(
  seurat_object = muscle,
  genes = genes_of_interest,
  grouping_column = "Injury", # Change to "Timepoint" or "Orig.ident" or as needed
  output_file = output_file_path
)

# Output: Access the pivot table for a specific gene if needed
View(gene_pivot_tables[["gene1"]])
