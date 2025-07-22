# Function: compare mNGS result species and label them with matched label (e.g. 1/n) .
# Author: Aixin LI (if you have issues, please contact liaixin@hku.hk)
# Update: 2025-06-02
# Description: matched label 1/n means: 1 for matched species between clinical and mNGS results, n for number of lab culture pathogens.


get_merged_short_table <- function(short_table, culture_table){
  # merge culture table with mNGS results
  short_table$Barcode_clean <- sprintf("%02d", as.integer(short_table$Barcode))
  culture_table <- culture_table %>%
    mutate(Barcode_clean = str_extract(str_remove_all(Barcode, "\\s+"), "\\d+$"))
  merged_short <- culture_table %>%
    left_join(short_table, by = c("Nanopore run", "Barcode_clean"), suffix = c("", "_mNGS")) %>%
    select(-Barcode_clean, -Barcode_mNGS)
  #merged_short[is.na(merged_short)] <- ""
  # **仅** 对字符列（character）把 NA 换成 ""
  merged_short <- merged_short %>%
    mutate(across(where(is.character), ~replace_na(.x, "")))

  return(merged_short)
}


get_merged_clean_table <- function(clean_table, culture_table, read_perct_cutoff=0){
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # merge culture table with mNGS results
  clean_table$Barcode <- sprintf("%02d", as.integer(clean_table$Barcode))
  
  report_list <- list()
  for (cutoff in read_perct_cutoff) {
    report <- clean_table %>%
      filter(Species != "") %>%
      filter(Adjusted_microbial_read_perc>=cutoff) %>%
      group_by(`Nanopore run`, Barcode, Database, Time, CS_score, Software) %>%
      arrange(desc(Adjusted_microbial_read_perc), .by_group = TRUE) %>%   
      mutate(n=n(), rank = paste0(row_number(), "/", n)) %>%
      summarise(mNGS = paste0(Species," (",Adjusted_microbial_read_perc,"%; ",rank, ")", collapse = "; ")) %>%
      ungroup() %>%
      mutate(Barcode_clean = str_extract(str_remove_all(Barcode, "\\s+"), "\\d+$")) %>%
      mutate(mNGS = gsub("; 1/1)",")", mNGS))
    
    report_list[[as.character(cutoff)]] <- report %>%
      mutate(CS_score = as.numeric(CS_score)) %>% 
      mutate(
        time_str = ifelse(is.na(Time), paste0(cutoff, "%"), paste0(Time, "hps_", cutoff, "%")),
        New_Column_Name = paste(
          paste0(Software, "-", Database), "Species",
          case_when(
            is.na(CS_score) ~ "",
            TRUE ~ paste0("cs", format(round(CS_score, 1), nsmall = 1))
          ),
          time_str, sep = "_"
        )
      ) %>%
      select(-CS_score,-Database,-Time,-Software,-time_str) %>%
      pivot_wider(names_from = New_Column_Name, values_from = mNGS, values_fn = list) %>%
      mutate(across(where(is.list), ~sapply(., function(x) paste(unlist(x), collapse = "; "))))
  }
  report <- reduce(report_list, full_join, by = c("Nanopore run", "Barcode", "Barcode_clean"))
  report[is.na(report)] <- "No species reported"
  
  culture_table <- culture_table %>%
    mutate(Barcode_clean = str_extract(str_remove_all(Barcode, "\\s+"), "\\d+$"))
  merged_clean <- culture_table %>%
    left_join(report, by = c("Nanopore run", "Barcode_clean"), suffix = c("", "_mNGS")) %>%
    select(-Barcode_clean, -Barcode_mNGS)
  merged_clean[is.na(merged_clean)] <- ""
  
  return(merged_clean)
}


add_match_label <- function(){
  # add match label
  # 参考BCB_Felice_v2.R的add match label，考虑隐藏功能
  
}
