#!/usr/bin/Rscript

# ------------------------------------------------------------
# Description: process Kraken-style reports (Kraken2, Bracken, etc.) into KT-format output
# ------------------------------------------------------------

# Author: Aixin LI
# Date: 30 May, 2025
# Note: If you have any questions, please contact liaixin@hku.hk


# 0. Load required packages ####
required_packages <- c("optparse", "stringr", "tools", "purrr", "dplyr", "tidyr", "tibble", "openxlsx", "readxl", "readr", "glue", "parallel")
invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  )
}))


# 1. Define command line options ####
option_list <- list(
  make_option(c("-f", "--filepath"), type = "character", default = NULL,
              help = "Input Nanopore folder, e.g. ./test/mixed_nanopore_run.", metavar = "character"),
  make_option(c("-c", "--culture_table"), type = "character", default = NULL, # Nanopore run, Barcode, Clinical culture result三列必需，且名称统一
              help = "[Optional] Clinical culture table in excel or tsv format.", metavar = "character"),
  make_option(c("--col_header_NR"), type = "character", default = NULL,  
              help = "[Optional] Specify the exact header of nanopore run column in culture table.", metavar = "character"),
  make_option(c("--col_header_Bc"), type = "character", default = NULL, 
              help = "[Optional] Specify the exact header of barcode column in culture table.", metavar = "character"),
  make_option(c("--col_header_Cl"), type = "character", default = NULL, 
              help = "[Optional] Specify the exact header of clinical culture result column in culture table.", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "[Optional] Output directory. Defaults: path/to/input/file/output_20XXXXXX-XXXX.", metavar = "character"),
  make_option(c("--read_count_cutoff"), type = "numeric", default = 1,
              help = "[Optional] Read count cutoff, input number: >=1. Default: 1.", metavar = "numeric"),
  make_option(c("--read_perct_cutoff"), type = "character", default = "0",
              help = "[Optional] Read percent cutoff (%), input number: 0-100. Comma-separated cutoff values, e.g. 10,20,30.", metavar = "cutoff_values"),
  make_option(c("--show_top"), type = "numeric", default = 100,
              help = "[Optional] Show top species/genus, input integer: >=1. Default: 100.", metavar = "numeric"),
  make_option(c("--time"), type = "character", default = "",
              help = "[Optional] Sequencing time (hps), e.g. 0.5/12/48/72 etc.", metavar = "character"),
  make_option(c("--threads"), type = "integer", default = 12,
            help = "[Optional] Number of threads to use for processing. Default: 12.")
)
# Parse command line arguments 
opt <- parse_args(OptionParser(option_list = option_list))


read_perct_cutoffs <- as.numeric(unlist(strsplit(opt$read_perct_cutoff, ",")))


if (!is.null(opt$filepath)) {
  if (is.null(opt$output)) {
    input_dir <- normalizePath(opt$filepath, mustWork = T)
    input_dir <- dirname(input_dir)
    
    opt$output <- paste0(input_dir,"/output_",format(Sys.time(), "%Y%m%d-%H%M"))
    dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
  }else{
    opt$output <- normalizePath(opt$output, mustWork = F)
  }
} else {
  stop("Error: input file (-f/--filepath) needed!")
}


start_time <- Sys.time()


# 2. Funciton defined ####
# script_dir <- dirname(normalizePath(sys.frame(1)$ofile)) # 只适用于Rstudio不适用Linux
get_script_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- grep("--file=", cmdArgs, value = TRUE)
  if (length(fileArg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", fileArg))))
  } else {
    # fallback: use current working directory
    return(getwd())
  }
}
script_dir <- get_script_dir()
source(file.path(script_dir, "Functions", "Handle_Kraken_Style_Reports.R"))
source(file.path(script_dir, "Functions", "Extract_Nanopore_Info.R"))
source(file.path(script_dir, "Functions", "Match_Judge.R"))


# 3. Run the function with parsed arguments ####
#设置folder批量读入，将结果汇总
message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Starting to read input Kraken2-style files")
message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Input file folder: ", opt$filepath)
file_list <- list.files(path = opt$filepath, pattern = "\\.kraken2_report$", full.names = TRUE)
num_cores <- opt$threads
processed_list <- mclapply(file_list, function(file_path) {
  file <- read_kraken_or_bracken(file_path = file_path)
  result <- handle_kraken_style_reports(
    df = file$data,
    software = file$type,
    
    read_count_cutoff = opt$read_count_cutoff,
    read_perct_cutoff = read_perct_cutoffs,
    top = opt$show_top,
    time = opt$time,
    
    run = extract_info(file_path)$run,
    barcode = extract_info(file_path)$barcode,
    database = extract_info(file_path)$database,
    CS_score = extract_info(file_path)$CS_score
  )
  
  # 调试: 给每个样本添加一个文件名列，后续检查
  sample_name <- basename(file_path)
  clean_df <- result$clean_table %>% mutate(filename = sample_name)
  short_df <- result$short_table %>% mutate(filename = sample_name)
  
  return(list(clean_df = clean_df, short_df = short_df))
}, mc.cores = num_cores)


message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Getting the clean and short mNGS result table")
merged_clean_df <- bind_rows(map(processed_list, "clean_df")) %>%
  mutate(Filename = filename) %>%
  select(-filename)
merged_short_df <- bind_rows(map(processed_list, "short_df")) %>%
  group_by(`Nanopore run`, Barcode) %>%
  summarise(across(everything(), ~ paste0(unique(.[!is.na(.) & . != ""]), collapse = "; "), .names = "{.col}"), .groups = "drop") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .))) %>%
  mutate(Filename = filename) %>%
  select(-filename)


# 4. Culture table merge (optional) ####
if (!is.null(opt$culture_table)) {
  message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Starting to merge mNGS result with culture table")
  message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Input culture result table: ", opt$culture_table)
  file_ext <- tools::file_ext(opt$culture_table)
  if (file_ext %in% c("xls", "xlsx", "xlsm")) {
    culture_table <- readxl::read_excel(opt$culture_table, sheet = 1)
  } else if (file_ext %in% c("tsv", "txt")) {
    culture_table <- readr::read_tsv(opt$culture_table)
  } else {
    stop("Error: unsupported file format for culture_table. Please provide a tab-separated (.tsv/.txt) or Excel (.xlsx/.xls) file.")
  }
  
  
  library(stringr)
  library(glue)
  
  rename_column <- function(df, pattern, new_name, user_defined_name = NULL, what = "parameter") {
    if (!is.null(user_defined_name)) {
      if (user_defined_name %in% names(df)) {
        names(df)[names(df) == user_defined_name] <- new_name
        message(glue::glue("Renamed user-specified column '{user_defined_name}' → '{new_name}'"))
      } else {
        stop(glue::glue("User-specified column '{user_defined_name}' not found in the table."))
      }
    } else {
      match_idx <- stringr::str_which(names(df), regex(pattern, ignore_case = TRUE))
      if (length(match_idx) == 1) {
        old_name <- names(df)[match_idx]
        names(df)[match_idx] <- new_name
        message(glue::glue("Renamed column '{old_name}' → '{new_name}' (matched by pattern '{pattern}')"))
      } else if (length(match_idx) == 0) {
        stop(glue::glue("No column '{new_name}' found. Please specify using --{what}"))
      } else {
        stop(glue::glue("Multiple columns match '{new_name}': {paste(names(df)[match_idx], collapse=', ')}. Please specify using --{what}"))
      }
    }
    return(df)
  }
  
  culture_table <- rename_column(culture_table, "nanopore[\\s_-]*run", "Nanopore run", user_defined_name = opt$col_header_NR, "col_header_NR")
  culture_table <- rename_column(culture_table, "barcode", "Barcode", user_defined_name = opt$col_header_Bc, "col_header_Bc")
  culture_table <- rename_column(culture_table, "clinical[\\s_-]*culture[\\s_-]*result", user_defined_name = opt$col_header_Cl, "Clinical culture result", "col_header_Cl")
  
  merged_result <- get_merged_short_table(merged_short_df, culture_table)
  # merged_result <- get_merged_clean_table(merged_clean_df, culture_table, read_perct_cutoffs)
}


# 5. Save outputs ####
input_name <- basename(opt$filepath)

if (!dir.exists(opt$output)) {
  message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Creating output dir:", opt$output)
  dir.create(opt$output, recursive = TRUE)
}

write.table(merged_clean_df, 
            file.path(opt$output, paste0("Table_Clean_mNGS_", input_name, "_", format(Sys.Date(), "%Y%m%d"), ".tsv")), 
            row.names = F, col.names = T, sep = "\t", quote = F)
write.table(merged_short_df, 
            file.path(opt$output, paste0("Table_Short_mNGS_", input_name, "_", format(Sys.Date(), "%Y%m%d"), ".tsv")), 
            row.names = F, col.names = T, sep = "\t", quote = F)

wb <- createWorkbook()
if(exists("merged_result")){
  addWorksheet(wb, "Merge result")
  writeData(wb, "Merge result", merged_result)
}
addWorksheet(wb, "Clean table")
writeData(wb, "Clean table", merged_clean_df)
addWorksheet(wb, "Short table")
writeData(wb, "Short table", merged_short_df)
saveWorkbook(wb, file.path(opt$output, paste0("Table_mNGS_", input_name, format(Sys.Date(), "%Y%m%d"), ".xlsx")), overwrite = TRUE)


end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Processing complete.")
message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Running time: ", round(elapsed, 4), " secs.")
message("[", format(Sys.time(), "%Y%m%d-%H:%M:%S"), "] Results saved in: ", opt$output, "\n")
