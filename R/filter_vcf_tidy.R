#' Filter VCF Data
#'
#' This function filters VCF data based on various criteria.
#'
#' @param input_file The path to the input file (.csv or .xlsx).
#' @param output_file The path to the output file (.csv).
#' @param zygosity_column The column name for zygosity filtering.
#' @param zygosity_value The value for zygosity filtering (e.g., "het" or "hom").
#' @param maf_criteria The criteria for MAF filtering (e.g., "<=0.001").
#' @param exome_column The column name for exomic allelic frequency in data. Missing values are replaced with 10.
#' @param genome_column The column name for genomic allelic frequency in data. Missing values are replaced with 10.
#' @param region_column The column name for region filtering.
#' @param region_filters A vector of filters for region filtering (e.g., c("exonic")).
#' @param effect_column The column name for effect filtering.
#' @param effect_filters A vector of filters for effect filtering (e.g., c("nonsynonymous SNV")).
#' @param gerp_column The column name for GERP filtering (optional).
#' @param gerp_criteria The criteria for GERP filtering (optional, e.g., ">=3").
#' @param cadd_column The column name for CADD filtering (optional).
#' @param cadd_criteria The criteria for CADD filtering (optional, e.g., ">=20").
#' @param custom_filters A list of custom filters (optional). Each filter is a list
#'   with two elements: `column` (the column name) and `criteria` (the filter criteria).
#' @return A filtered data frame.
#' @examples
#'\dontrun{
#' filter_vcf_tidy(
#'   input_file = "your_input_file.csv",
#'   output_file = "your_output_file.csv",
#'   zygosity_column = "Zygosity",
#'   zygosity_value = "het",
#'   maf_criteria = "<=0.001",
#'   exome_column = "gnomAD_exome_ALL",
#'   genome_column = "gnomAD_genome_ALL",
#'   region_column = "Func.refGene",
#'   region_filters = c("exonic"),
#'   effect_column = "ExonicFunc.refGene",
#'   effect_filters = c("nonsynonymous SNV", "stoploss"),
#'   gerp_column = "GERP.._RS",
#'   gerp_criteria = ">=3",
#'   custom_filters = list(
#'     list(column = "SIFT_pred", criteria = "== 'D'"),
#'     list(column = "PROVEAN_pred", criteria = "== 'D'"),
#'     list(column = "LRT_pred", criteria = "== 'D'")
#'   )
#' )
#' }
#' @importFrom dplyr %>% mutate across everything rename_with all_of if_else filter .data
#' @importFrom vroom vroom cols
#' @importFrom readxl read_excel
#' @importFrom readr write_csv
#' @export
filter_vcf_tidy <- function(input_file, output_file, zygosity_column, zygosity_value, maf_criteria,
                            exome_column, genome_column, region_column, region_filters,
                            effect_column, effect_filters, gerp_column = NULL, gerp_criteria = NULL,
                            cadd_column = NULL, cadd_criteria = NULL, custom_filters = NULL) {

  # Read the input file based on its extension
  file_extension <- tools::file_ext(input_file)
  if (file_extension == "csv") {
    vcf_data <- vroom(input_file, col_types = cols(.default = "c"))
  } else if (file_extension == "xlsx") {
    vcf_data <- read_excel(input_file) %>% mutate(across(everything(), as.character))
  } else {
    stop("Unsupported file format. Please provide a .csv or .xlsx file.")
  }

  # Ensure unique and valid column names
  vcf_data <- vcf_data %>% rename_with(~make.names(., unique = TRUE))

  # Print column names to verify
  print("Column names in the dataset:")
  print(colnames(vcf_data))

  # Check for required columns
  required_columns <- c(region_column, effect_column, zygosity_column, exome_column, genome_column)
  if (!is.null(gerp_column)) {
    required_columns <- c(required_columns, gerp_column)
  }
  if (!is.null(cadd_column)) {
    required_columns <- c(required_columns, cadd_column)
  }

  missing_columns <- setdiff(required_columns, colnames(vcf_data))
  if (length(missing_columns) > 0) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  # Safely convert columns to numeric and replace NAs
  safe_numeric_conversion <- function(column, default_na_value) {
    as.numeric(dplyr::if_else(is.na(column) | column == "", as.character(default_na_value), as.character(column)))
  }

  vcf_data <- vcf_data %>%
    mutate(across(all_of(exome_column), ~safe_numeric_conversion(., 10), .names = "clean_{col}")) %>%
    mutate(across(all_of(genome_column), ~safe_numeric_conversion(., 10), .names = "clean_{col}")) %>%
    mutate(across(all_of(gerp_column), ~safe_numeric_conversion(., -1), .names = "clean_{col}")) %>%
    mutate(across(all_of(cadd_column), ~safe_numeric_conversion(., -40), .names = "clean_{col}"))

  # Filtering steps
  # Region filtering
  vcf_data <- vcf_data %>%
    filter(.data[[region_column]] %in% region_filters)
  print("Data summary after region filtering:")
  print(dim(vcf_data))

  # Effect filtering
  vcf_data <- vcf_data %>%
    filter(.data[[effect_column]] %in% effect_filters)
  print("Data summary after effect filtering:")
  print(dim(vcf_data))

  # Zygosity filtering
  vcf_data <- vcf_data %>%
    filter(.data[[zygosity_column]] == zygosity_value)
  print("Data summary after zygosity filtering:")
  print(dim(vcf_data))

  # Convert exome and genome columns to numeric and filter based on MAF criteria
  exome_operator <- substr(maf_criteria, 1, 2)
  exome_threshold <- as.numeric(substr(maf_criteria, 3, nchar(maf_criteria)))
  genome_operator <- substr(maf_criteria, 1, 2)
  genome_threshold <- as.numeric(substr(maf_criteria, 3, nchar(maf_criteria)))

  vcf_data <- vcf_data %>%
    filter(eval(parse(text = paste0("clean_", exome_column, exome_operator, exome_threshold)))) %>%
    filter(eval(parse(text = paste0("clean_", genome_column, genome_operator, genome_threshold))))
  print("Data summary after MAF filtering:")
  print(dim(vcf_data))

  # Filter based on GERP criteria
  if (!is.null(gerp_column) && !is.null(gerp_criteria)) {
    gerp_operator <- substr(gerp_criteria, 1, 2)
    gerp_threshold <- as.numeric(substr(gerp_criteria, 3, nchar(gerp_criteria)))
    vcf_data <- vcf_data %>%
      filter(eval(parse(text = paste0("clean_", gerp_column, gerp_operator, gerp_threshold))))
    print("Data summary after GERP filtering:")
    print(dim(vcf_data))
  }

  # Filter based on CADD criteria
  if (!is.null(cadd_column) && !is.null(cadd_criteria)) {
    cadd_operator <- substr(cadd_criteria, 1, 2)
    cadd_threshold <- as.numeric(substr(cadd_criteria, 3, nchar(cadd_criteria)))
    vcf_data %>%
      filter(eval(parse(text = paste0("clean_", cadd_column, cadd_operator, cadd_threshold))))
    print("Data summary after CADD filtering:")
    print(dim(vcf_data))
  }

  # Custom filters
  if (!is.null(custom_filters)) {
    for (filter in custom_filters) {
      custom_column <- filter$column
      custom_criteria <- filter$criteria
      vcf_data <- vcf_data %>%
        filter(eval(parse(text = paste0(".data[['", custom_column, "']] ", custom_criteria))))
      print(paste("Data summary after custom filtering on column:", custom_column))
      print(dim(vcf_data))
    }
  }

  # Write the filtered data to a CSV file
  write_csv(vcf_data, output_file)
  message("Filtered data has been written to ", output_file)
}
