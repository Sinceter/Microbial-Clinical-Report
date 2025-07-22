# Function: read kraken-style format and rearrange them into readable format
# Author: Aixin LI (if you have issues, please contact liaixin@hku.hk)
# Update: 2025-05-30
# Description: support Kraken2 PlusPF database version 2025-04-02, in this version, Viruses belong to Domain R1 whereas Bacteria/Archaea/Eukaryota belong to Domain R2. However, in version 2023-03-14, they all belong to Domain D. Therefore, we use taxonomy ID to locate the four domains (Viruses/Bacteria/Archaea/Eukaryota).
# 用于merge culture table的：kraken2$short_table
# 用于中间处理的：kraken2$clean_table
## No microbial reads -> not microbial reads found
## No species reported -> microbial reads found, whereas no species can be reported under this threshold combination


read_kraken_or_bracken <- function(file_path) {
  kraken2_header <- c("% of Seqs", "Clades", "Taxonomies", "Kmers", "Distinct Kmers", "Rank", "Taxonomy ID", "Scientific Name")
  bracken_header <- c("% of Seqs", "Clades", "Taxonomies", "Rank", "Taxonomy ID", "Scientific Name")
  
  first_line <- readLines(file_path, n = 1)
  split_fields <- unlist(strsplit(first_line, "\t"))
  
  if (all(kraken2_header %in% split_fields)) {
    df <- read.delim(file_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    file_type <- "Kraken2"
  } else if (all(bracken_header %in% split_fields)) {
    df <- read.delim(file_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    file_type <- "Bracken"
    df$Kmers <- NA
    df$`Distinct Kmers` <- NA
  } else {
    # 没有表头，强制读取为 Bracken 格式
    df <- read.delim(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(df) == 6) {
      colnames(df) <- bracken_header
      df$Kmers <- NA
      df$`Distinct Kmers` <- NA
      file_type <- "Bracken"
    } else if (ncol(df) == 8) {
      colnames(df) <- kraken2_header
      file_type <- "Kraken2"
    } else {
      stop("Unknown file format: cannot determine column names.")
    }
  }
  
  return(list(data = df, type = file_type))
}


handle_kraken_style_reports <- function(df, 
                                    read_count_cutoff=1, read_perct_cutoff=0, top=100,
                                    software = "", # Kraken2 or Bracken
                                    run = "",
                                    time = "",
                                    barcode = "",
                                    database="",
                                    CS_score = "",
                                    rm_HM=T){ 

  # re-organize the input file
  lines <- split(df, seq(nrow(df)))  # 每行作为列表
  
  
  tmp <- data.frame(`Nanopore run` = character(), Time = character(),
                    Barcode = character(), Database = character(),
                    CS_score = character(), Software = character(), 
                    Domain = character(), Kingdom = character(), Phylum = character(),
                    Class = character(), Order = character(), Family = character(),
                    Genus = character(), subGenus = character(), Species = character(),
                    Read_counts = numeric(), Read_percent_from_report = numeric(),
                    stringsAsFactors = FALSE, check.names = FALSE)
  
  
  valid_ranks <- df$Rank[df$Rank %in% c("K", "P", "C", "O", "F", "G", "G1", "S")]
  if (length(valid_ranks) == 0) { # 处理100%为unassigned的kraken2 report
    # return empty structure
    empty_clean <- tmp
    empty_short <- data.frame(`Nanopore run` = run, Barcode = barcode,
                              stringsAsFactors = FALSE, check.names = FALSE)
    # 短表(KT表)表头
    for (cutoff in read_perct_cutoff) {
      # 构建动态列名
      time_suffix <- if (time == "") paste0(cutoff, "%") else paste0(time, "hps_", cutoff, "%")
      gheader <- paste(if (database == "") software else paste0(software, "-", database),
                       "Genus",
                       if (CS_score != "") paste0("cs", CS_score) else NULL,
                       time_suffix,
                       sep = "_")
      sheader <- paste(if (database == "") software else paste0(software, "-", database),
                       "Species",
                       if (CS_score != "") paste0("cs", CS_score) else NULL,
                       time_suffix,
                       sep = "_")
      
      empty_short[[gheader]] <- "No microbial reads"
      empty_short[[sheader]] <- "No microbial reads"
    }

    return(list(clean_table = empty_clean, short_table = empty_short)) # 立即结束，剩余代码不会被运行
  }
  
  
  # init the indices
  read_count <- ""
  level_indx <- c("domain","kingdom","phylum","class","order","family","genus","subgenus")
  for (name in level_indx) {
    assign(name, "")  
  }
  num_indx <- c("d","k","p","c","o","f","g","g1","s") # 9
  for (n in num_indx) {
    assign(n, NA)  
  }
  
  # process the file row by row
  for (line in lines) {
    # line <- trimws(line)  # trim starting spaces
    # fields <- strsplit(line, "\t")[[1]]
    # read_perc <- fields[1]   # % of Seqs
    # read_count <- fields[2]  # Clades
    # rank <- fields[6] # Rank
    # name <- trimws(gsub("^[ ]*", "", fields[8])) # Scientific Name
    # leading_spaces <- sapply(fields[8], function(x) nchar(str_extract(x, "^\\s*")))
    library(stringr)
    read_perc <- as.character(line[["% of Seqs"]])
    read_count <- as.character(line[["Clades"]])
    rank <- as.character(line[["Rank"]])
    taxonID <- as.character(line[["Taxonomy ID"]]) # 
    raw_name <- as.character(line[["Scientific Name"]])
    name <- trimws(gsub("^[ ]*", "", raw_name))
    leading_spaces <- sapply(raw_name, function(x) nchar(str_extract(x, "^\\s*")))
    
    
    if (taxonID %in% c("10239", "2", "2157", "2759")) { # they are: Viruses/Bacteria/Archaea/Eukaryota
      domain <- name
      d <- leading_spaces
      for (name in level_indx[-1]) { assign(name, "") }
      for (n in num_indx[-1]) { assign(n, NA) }
    } else if (rank == "R") {  # Root_minus_human reads, based on all classified reads
      kraken_root_reads <- as.numeric(read_count)
    } else if (rank == "K") {
      kingdom <- name
      k <- leading_spaces
      if (!is.na(d) & k<=d) { domain <- "" }
      for (name in level_indx[-(1:2)]) { assign(name, "") }
      for (n in num_indx[-(1:2)]) { assign(n, NA) }
    } else if (rank == "P") {
      phylum <- name
      p <- leading_spaces # p>=k+2
      if (!is.na(k) & p<=k) { 
        kingdom <- ""
        if (!is.na(d) & p<=d) { domain <- "" }
      }
      for (name in level_indx[-(1:3)]) { assign(name, "") }
      for (n in num_indx[-(1:3)]) { assign(n, NA) }
    } else if (rank == "C") {
      class <- name
      c <- leading_spaces 
      if (!is.na(p) & c<=p){
        phylum <- ""
        if (!is.na(k) & c<=k) { 
          kingdom <- ""
          if (!is.na(d) & c<=d) { domain <- "" }
        }
      }
      for (name in level_indx[-(1:4)]) { assign(name, "") }
      for (n in num_indx[-(1:4)]) { assign(n, NA) }
    } else if (rank == "O") {
      order <- name
      o <- leading_spaces 
      if (!is.na(c) & o<=c) {
        class <- ""
        if (!is.na(p) & o<=p){
          phylum <- ""
          if (!is.na(k) & o<=k) { 
            kingdom <- ""
            if (!is.na(d) & o<=d) { domain <- "" }
          }
        }
      }
      for (name in level_indx[-(1:5)]) { assign(name, "") }
      for (n in num_indx[-(1:5)]) { assign(n, NA) }
    } else if (rank == "F") {
      family <- name
      f <- leading_spaces 
      if (!is.na(o) & f<=o) {
        order <- ""
        if (!is.na(c) & f<=c) {
          class <- ""
          if (!is.na(p) & f<=p){
            phylum <- ""
            if (!is.na(k) & f<=k) { 
              kingdom <- ""
              if (!is.na(d) & f<=d) { domain <- "" }
            }
          }
        }
      }
      for (name in level_indx[-(1:6)]) { assign(name, "") }
      for (n in num_indx[-(1:6)]) { assign(n, NA) }
    } else if (rank == "G") {
      genus <- name
      g <- leading_spaces 
      if (!is.na(f) & g<=f) {
        family <- ""
        if (!is.na(o) & g<=o) {
          order <- ""
          if (!is.na(c) & g<=c) {
            class <- ""
            if (!is.na(p) & g<=p){
              phylum <- ""
              if (!is.na(k) & g<=k) { 
                kingdom <- ""
                if (!is.na(d) & g<=d) { domain <- "" }
              }
            }
          }
        }
      }
      for (name in level_indx[-(1:7)]) { assign(name, "") }
      for (n in num_indx[-(1:7)]) { assign(n, NA) }
      tmp <- rbind(tmp, data.frame(`Nanopore run` = run, Time = time,
                                   Barcode = barcode, Database = database,
                                   CS_score = CS_score, Software = software, 
                                   Domain = domain, Kingdom = kingdom, Phylum = phylum,
                                   Class = class, Order = order, Family = family,
                                   Genus = genus, subGenus = "", Species = "", 
                                   Read_counts = as.numeric(read_count), 
                                   Read_percent_from_report = as.numeric(read_perc), # Read_percent_from_report = numeric(), # Root_minus_human reads, based on all classified reads
                                   stringsAsFactors = FALSE, check.names = FALSE))
    } else if (rank == "G1") {
      subgenus <- name
      g1 <- leading_spaces 
      if (!is.na(g) & g1<=g) {
        genus <- ""
        if (!is.na(f) & g1<=f) {
          family <- ""
          if (!is.na(o) & g1<=o) {
            order <- ""
            if (!is.na(c) & g1<=c) {
              class <- ""
              if (!is.na(p) & g1<=p){
                phylum <- ""
                if (!is.na(k) & g1<=k) { 
                  kingdom <- ""
                  if (!is.na(d) & g1<=d) { domain <- "" }
                }
              }
            }
          }
        }
      }
      for (name in level_indx[-(1:8)]) { assign(name, "") }
      for (n in num_indx[-(1:8)]) { assign(n, NA) }
      tmp <- rbind(tmp, data.frame(`Nanopore run` = run, Time = time,
                                   Barcode = barcode, Database = database,
                                   CS_score = CS_score, Software = software, 
                                   Domain = domain, Kingdom = kingdom, Phylum = phylum,
                                   Class = class, Order = order, Family = family,
                                   Genus = genus, subGenus = subgenus, Species = "", 
                                   Read_counts = as.numeric(read_count), 
                                   Read_percent_from_report = as.numeric(read_perc),
                                   stringsAsFactors = FALSE, check.names = FALSE))
    } else if (rank == "S") {
      species <- name
      s <- leading_spaces 
      if (!is.na(g1) & s<=g1) {
        subgenus <- ""
        if (!is.na(g) & s<=g) {
          if (!is.na(f) & s<=f) {
            family <- ""
            if (!is.na(o) & s<=o) {
              order <- ""
              if (!is.na(c) & s<=c) {
                class <- ""
                if (!is.na(p) & s<=p){
                  phylum <- ""
                  if (!is.na(k) & s<=k) { 
                    kingdom <- ""
                    if (!is.na(d) & s<=d) { domain <- "" }
                  }
                }
              }
            }
          }
        }
      }
      tmp <- rbind(tmp, data.frame(`Nanopore run` = run, Time = time,
                                   Barcode = barcode, Database = database,
                                   CS_score = CS_score, Software = software, 
                                   Domain = domain, Kingdom = kingdom, Phylum = phylum,
                                   Class = class, Order = order, Family = family,
                                   Genus = genus, subGenus = subgenus, Species = species, 
                                   Read_counts = as.numeric(read_count), 
                                   Read_percent_from_report = as.numeric(read_perc),
                                   stringsAsFactors = FALSE, check.names = FALSE))
    }
  }
  
  # Root_minus_human reads, based on all classified reads
  human_reads <- tmp[which(tmp$Species=="Homo sapiens"),"Read_counts"]
  human_reads <- ifelse(length(human_reads)==0, 0, human_reads)
  nonhuman_root_reads <- as.numeric(kraken_root_reads)-as.numeric(human_reads)
  tmp$Adjusted_microbial_read_perc <- round((tmp$Read_counts/nonhuman_root_reads)*100,2) # read percentage based on all microbial reads
  
  # remove Homo sapiens
  filter_tmp <- tmp
  if (rm_HM) {
    filter_tmp <- tmp[which(tmp$Genus!="Homo"),]
  }
  
  # criteria to be reported: >= read_count_cutoff & >= read_perc_cutoff
  filter_tmp <- filter_tmp[which(filter_tmp$Read_counts >= read_count_cutoff),]
  # filter_tmp <- filter_tmp[which(filter_tmp$Read_percent_from_report >= read_perct_cutoff),]
  
  # KT format
  filter_tmp <- filter_tmp[order(-filter_tmp$Adjusted_microbial_read_perc),]
  # species_only <- filter_tmp[filter_tmp$Species!="",]
  # genus_only <- filter_tmp[filter_tmp$subGenus=="" & filter_tmp$Species=="",
  #                          -c(which(colnames(filter_tmp)=="subGenus"), which(colnames(filter_tmp)=="Species"))]

  # top_genus <- paste(paste0(genus_only$Genus[1:top], " (", genus_only$Adjusted_microbial_read_perc[1:top], "%, ", genus_only$Read_counts[1:top],")"), collapse = "; ")
  # top_genus <- gsub(";?\\s*NA \\(NA%, NA\\)", "", top_genus)
  # top_genus <- ifelse(top_genus=="", "No species reported", top_genus)
  # top_spp <- paste(paste0(species_only$Species[1:top], " (", species_only$Adjusted_microbial_read_perc[1:top], "%, ", species_only$Read_counts[1:top],")"), collapse = "; ")
  # top_spp <- gsub(";?\\s*NA \\(NA%, NA\\)", "", top_spp)
  # top_spp <- ifelse(top_spp=="", "No species reported", top_spp)
  
  short_table <- data.frame(
    `Nanopore run` = run, Barcode = barcode,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  for (cutoff in read_perct_cutoff) {
    filter_cut <- filter_tmp[filter_tmp$Adjusted_microbial_read_perc >= cutoff, ]
    species_only <- filter_cut[filter_cut$Species != "", ]
    genus_only <- filter_cut[filter_cut$subGenus == "" & filter_cut$Species == "", 
                             -c(which(colnames(filter_cut) == "subGenus"), which(colnames(filter_cut) == "Species"))]
    
    top_genus <- paste(paste0(genus_only$Genus[1:top], " (", genus_only$Adjusted_microbial_read_perc[1:top], "%, ", genus_only$Read_counts[1:top], ")"), collapse = "; ")
    top_genus <- gsub(";?\\s*NA \\(NA%, NA\\)", "", top_genus)
    top_genus <- ifelse(top_genus == "", "No genus reported", top_genus)
    top_spp <- paste(paste0(species_only$Species[1:top], " (", species_only$Adjusted_microbial_read_perc[1:top], "%, ", species_only$Read_counts[1:top], ")"), collapse = "; ")
    top_spp <- gsub(";?\\s*NA \\(NA%, NA\\)", "", top_spp)
    top_spp <- ifelse(top_spp == "", "No species reported", top_spp)
    
    time_suffix <- if (time == "") paste0(cutoff, "%") else paste0(time, "hps_", cutoff, "%")
    gheader <- paste(if (database == "") software else paste0(software, "-", database),
                     "Genus",
                     if (CS_score != "") paste0("cs", CS_score) else NULL,
                     time_suffix,
                     sep = "_")
    sheader <- paste(if (database == "") software else paste0(software, "-", database),
                     "Species",
                     if (CS_score != "") paste0("cs", CS_score) else NULL,
                     time_suffix,
                     sep = "_")
    
    short_table[[gheader]] <- top_genus
    short_table[[sheader]] <- top_spp
  }

  return(list(clean_table=filter_tmp,
              short_table=short_table))
}
