# Function: extract nanopore/barcode/database/cs score information from set filename
# Author: Aixin LI (if you have issues, please contact liaixin@hku.hk)
# Update: 2025-05-31
# Description: to make it more convenient for user skipping the step of manually inputting above-mentioned information.


extract_info <- function(file_path) {
  filename <- basename(file_path)
  
  run <- paste0("Nanopore", sub("_nohuman.*", "", filename)) # 默认是没有Nanopore几个字母的

  barcode <- sub(".*barcode(\\d+).*", "\\1", filename)

  CS_score <- sub(".*\\.(\\d+\\.\\d+).*", "\\1", filename)
  
  # 匹配 .<db>.<cs_score>
  database <- sub(paste0(".*\\.([^.]+)\\.", CS_score, ".*"), "\\1", filename)
  
  return(list(
    path = file_path,
    run = run,
    barcode = barcode,
    database = database,
    CS_score = CS_score
  ))
}
