test_that("filter_vcf_tidy handles missing values correctly", {
  # Create a temporary CSV file with example data
  temp_file <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    Zygosity = c("het", "hom", "het", "hom"),
    gnomAD_exome_ALL = c(NA, 0.002, 0.005, NA),
    gnomAD_genome_ALL = c(0.001, NA, 0.003, 0.004),
    Func.refGene = c("exonic", "intronic", "exonic", "exonic"),
    ExonicFunc.refGene = c("nonsynonymous SNV", "synonymous SNV", "nonsynonymous SNV", "stoploss"),
    GERP.._RS = c(4, 2, 3.5, NA),
    SIFT_pred = c("D", "T", "D", "D")
  ), temp_file, row.names = FALSE)

  output_file <- tempfile(fileext = ".csv")

  filter_vcf_tidy(
    input_file = temp_file,
    output_file = output_file,
    zygosity_column = "Zygosity",
    zygosity_value = "het",
    maf_criteria = "<=0.001",
    exome_column = "gnomAD_exome_ALL",
    genome_column = "gnomAD_genome_ALL",
    region_column = "Func.refGene",
    region_filters = c("exonic"),
    effect_column = "ExonicFunc.refGene",
    effect_filters = c("nonsynonymous SNV", "stoploss"),
    gerp_column = "GERP.._RS",
    gerp_criteria = ">=3",
    custom_filters = list(
      list(column = "SIFT_pred", criteria = "== 'D'")
    )
  )

  result <- read.csv(output_file)

  # Check if missing values were replaced with 10
  expect_true(all(result$clean_gnomAD_exome_ALL >= 10 | !is.na(result$clean_gnomAD_exome_ALL)))
  expect_true(all(result$clean_gnomAD_genome_ALL >= 10 | !is.na(result$clean_gnomAD_genome_ALL)))
})

