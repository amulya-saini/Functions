# Load the required libraries
library(biomaRt)
library(writexl)
library(readxl)

#' Convert Gene Names Between Species Using Ensembl Biomart
#'
#' This function reads a gene dataset from an Excel file, 
#' retrieves orthologous gene names for a target species using Ensembl Biomart, 
#' and writes the converted data to a new Excel file or a new sheet in the same file.
#'
#' @param input.file Character. Path to the input Excel file containing gene names.
#' @param sheet.name Character. Name of the sheet in the Excel file to process.
#' @param species.from Character. Source species in Biomart format (e.g., "mmusculus" for mouse).
#' @param species.to Character. Target species in Biomart format (e.g., "hsapiens" for human).
#' @param output.file Character. Path to save the output Excel file with converted gene names.
#' @param new.file Logical. If TRUE, writes to a new file. If FALSE, adds a new sheet to the input file.
#' @param new.sheet.name Character. Name of the new sheet to add if new.file is FALSE.
#' 
#' @input.format The input data should have cell types as column names 
#' and markers/gene names as row values.
#'
#' Example format:
#' 
#' | Cell_Type1 | Cell_Type2 |
#' |------------|------------|
#' | GeneA      | GeneX      |
#' | GeneB      | GeneY      |
#' | GeneC      |            |
#' |            | GeneZ      |
#'
#' @return None. The function writes the output to a specified Excel file or appends it as a new sheet.
#'
#' @examples
#' convert_species_genes("Zebrafish_Cell_type_markers_db.xlsx", 
#'                      "Immune.cell.zebrafish", 
#'                      "mmusculus", "hsapiens", 
#'                      "Human.immune.Genes.Mapped.from.Immune.zebrafish.xlsx", 
#'                      FALSE)
#'
#' convert_species_genes("Zebrafish_Cell_type_markers_db.xlsx", 
#'                      "Immune.cell.zebrafish", 
#'                      "mmusculus", "hsapiens",
#'                      new.file = TRUE, new.sheet.name = "Converted_Genes")                       
#'                      
convert_species_genes <- function(input.file, 
                                  sheet.name, 
                                  species.from, 
                                  species.to, 
                                  output.file, 
                                  new.file = TRUE, 
                                  new.sheet.name = "Converted_Genes") {
  
  # Load the Excel file and the relevant sheet
  data <- read_excel(input.file, sheet = sheet.name)
  
  # Establish connections to Ensembl Biomart
  ensembl_from <- useEnsembl(biomart = "ensembl", 
                             dataset = paste0(species.from, "_gene_ensembl"), 
                             version = 112) # Update the version number based on recent ensembl archive
                                            # https://useast.ensembl.org/info/website/archives/index.html
  
  # Function to retrieve orthologues for the specified species
  get_orthologues <- function(genes) {
    
    orthologues <- getBM(
      attributes = c("external_gene_name", 
                     paste0(species.to, "_homolog_associated_gene_name")),
      filters = "external_gene_name",
      values = genes,
      mart = ensembl_from,
      uniqueRows = TRUE)
    
    return(orthologues)
  }
  
  # Process each column in the data
  converted_data <- data
  for (col in colnames(data)) {
    original_genes <- na.omit(data[[col]]) # Remove NA values
    orthologues <- get_orthologues(original_genes)
    converted_genes <- orthologues[[paste0(species.to, "_homolog_associated_gene_name")]]
    
    # Replace original genes with corresponding converted genes
    converted_data[[col]] <- ifelse(
      data[[col]] %in% orthologues$external_gene_name,
      orthologues[[paste0(species.to, "_homolog_associated_gene_name")]][match(data[[col]], orthologues$external_gene_name)],NA)
  }
  
  # Save the result
  if (new.file) {
    if (is.null(output.file)) {
      stop("Please provide a valid output.file path when new.file = TRUE.")
    }
    # Write to a new file
    write_xlsx(converted_data, output.file)
  } else {
    # Append as a new sheet to the input file
    existing_sheets <- excel_sheets(input.file)
    all_sheets <- lapply(existing_sheets, function(sheet) read_excel(input.file, sheet = sheet))
    names(all_sheets) <- existing_sheets
    all_sheets[[new.sheet.name]] <- converted_data
    
    # Overwrite the input file with the new sheet added
    write_xlsx(all_sheets, input.file)
  }
}
