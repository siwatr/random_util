#!/usr/bin/env Rscript

# Deeptools command generator:
# Generate command for DeepTools' computeMatrix and run it.

#' format_argv
#'
#' This function processes a vector of command-line arguments (argv) and organizes them into a list format.
#' It identifies flags (arguments starting with '-') and classifies them based on the number of leading dashes.
#' The function also counts the number of arguments associated with each flag and can optionally keep or discard
#' leading non-flag arguments.
#'
#' @param argv A character vector of command-line arguments.
#' @param keep_leading_info A logical value indicating whether to keep leading non-flag arguments in the output. Default is FALSE.
#'
#' @return A list containing two elements:
#' \item{info}{A data frame with details about each argument, including its name, flag type, and the number of associated arguments.}
#' \item{arg}{A list where each element corresponds to a flag and contains its associated arguments.}
#'
#' @import dplyr
#' @importFrom stringr str_detect str_remove str_count str_extract
#' @examples
#' argv <- c("input.txt", "-a", "value1", "--long-flag", "value2", "-b")
#' result <- format_argv(argv)
#' print(result)
#'

format_argv <- function(argv, keep_leading_info = FALSE, keep_dash = FALSE) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
  })
  # Turn argv vector into a list
  arg_name_flag <- str_detect(argv, "^\\-")

  args_info <- data.frame(received_arg = c("NA", argv[arg_name_flag]))
  # temporary function to classify flag type
  classify_flag_type <- function(dash_count) {
    # dash_count <- c(NA, 0, 1, 2, 2, 1, 3, 4, 5 ,1)
    if (any(is.na(dash_count))) {
      dash_count[is.na(dash_count)] <- 0
    }

    dash_count <- pmin(dash_count, 3)

    flag_type <- sapply(dash_count, FUN = function(x) {
      switch(paste0("f", x),
        f0 = "leading", # First, non-flag argument
        f1 = "short",
        f2 = "long",
        f3 = "unknown"
      )
    })
    return(flag_type)
  }

  args_info <- args_info %>%
    mutate(
      arg_name = str_remove(received_arg, "^\\-+"),
      flag_dash_count = str_count(str_extract(received_arg, "^\\-+"), "\\-"),
      flag_type = classify_flag_type(flag_dash_count),
      is_flag = FALSE,
      n_args = 0
    )

  # Organize cmd argument into list object
  args_list <- list()
  cur_arg_name <- NULL
  n_arg <- 0
  for (i in seq_along(argv)) {
    if (arg_name_flag[i]) {
      n_arg <- n_arg + 1
      cur_arg_name <- str_remove(argv[i], "^\\-+")

      is_flag <- FALSE
      if (i == length(argv)) {
        is_flag <- TRUE
      } else if (arg_name_flag[i + 1]) {
        is_flag <- TRUE
      }

      if (is_flag) {
        args_list[[cur_arg_name]] <- TRUE
        args_info$is_flag[args_info$arg_name == cur_arg_name] <- TRUE
        args_info$n_args[args_info$arg_name == cur_arg_name] <- 1
      }
    } else {
      if (is.null(cur_arg_name)) {
        # i.e., leading arg
        if (length(args_list) == 0) {
          args_list[[1]] <- argv[i]
        } else {
          args_list[[1]] <- c(args_list[[1]], argv[i])
        }
        args_info$n_args[1] <- args_info$n_args[1] + 1
      } else {
        # add value to flag
        args_list[[cur_arg_name]] <- c(args_list[[cur_arg_name]], argv[i])
        args_info$n_args[args_info$arg_name == cur_arg_name] <- args_info$n_args[args_info$arg_name == cur_arg_name] + 1
      }
    }
  }

  if (!keep_leading_info) {
    if (args_info$n_args[1] == 0) {
      # i.e., leading args
      args_info[-1, ]
    }
  }

  if(keep_dash){
    names(args_list) <- replace(args_info$received_arg, is.na(args_info$received_arg), "leading_args")
  }

  return(list(
    info = args_info,
    arg = args_list
  ))
}

split_args_txt <- function(args_txt, keep_dash=FALSE) {
  # mainly use this one for debug
  # args_txt <- " --asdf x y z -a 1 2 -z -xy --XY -b 3 4 5 "
  require(stringr)
  args_txt <- str_remove(args_txt, "(^\\s+)|(\\s$)")
  args_vec <- str_split_1(args_txt, "\\s+(?=\\-)")

  # Identify flag type
  dash_count <- str_count(args_vec, "(^\\-)|((?<=\\-)\\-)")
  args_type <- sapply(dash_count, function(x) {
    if (x == 0) {
      return("leading")
    } else if (x == 1) {
      return("short")
    } else if (x == 2) {
      return("long")
    } else {
      return("unknown")
    }
  })

  # Unpack args and add to args_list
  args_list <- list()
  for (i in seq_along(args_vec)) {
    A = str_split_1(args_vec[i], "\\s+")
    A = A[A != ""] # just in case
    if (args_type[i] == "leading") {
      args_list == c(args_list, list(leading = A))
    } else {
      # named args
      n_args = length(A) - 1
      if (n_args == 0) {
        # i.e., a flag
        if (args_type[i] == "short") {
          f <- str_split_1(A, "")
          f <- f[f != "-"]
          if(keep_dash) f <- paste0("-", f)
          f_list <- as.list(rep(TRUE, length(f)))
          names(f_list) <- f
          args_list <- c(args_list, f_list)
        } else {
          A_list <- list(TRUE)
          names(A_list) <- ifelse(keep_dash, yes = A, no = str_remove(A, "^\\-+"))
          args_list <- c(args_list, A_list)
        }
      } else {
        # i.e., not a flag
        A_list <- list(A[-1])
        names(A_list) <- ifelse(keep_dash, yes = A[1], no = str_remove(A[1], "^\\-+"))
        args_list <- c(args_list, A_list)
      }
    } # end check leading
  } # end loop

  return(args_list)
}
# Quick test
if(F) split_args_txt(" --asdf x y z -a 1 2 -z -xy --XY -b 3 4 5 ")

# Format sample (e.g., BAM or bigwig) inputs ======================================================
#' Convert a table to command-line arguments
#'
#' This function reads a table from a specified file path and converts its contents into a vector of command-line arguments.
#' It allows for optional filtering and reordering of entries, as well as renaming of argument names.
#'
#' @param table_path A string specifying the path to the input table file.
#' @param use_entry A vector of entries (i.e., sample) to use from the table. If NULL, all entries from the first column are used.
#' @param use_args A named vector of arguments to use. If NULL, all arguments are used.
#' @param collapse_redundant A logical value indicating whether to collapse redundant argument values. Default is TRUE.
#' @return A character vector of command-line arguments.
#' @examples
#' table_to_args("path/to/table.txt", use_entry = c("entry1", "entry2"), use_args = c("-arg1", "-arg2"))
#'
table_to_args <- function(table_path, use_entry = NULL, use_args = NULL, collapse_redundant = TRUE, as_text=FALSE){
  ref_df <- read.table(table_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if(is.null(use_entry)){
    # Which rows of the table to be used
    use_entry <- ref_df[[1]]
  }

  # ref_df <- dplyr::filter(ref_df, label %in% use_entry)
  ref_df <- ref_df[ref_df[[1]] %in% use_entry, ]
  ref_df <- ref_df[match(use_entry, ref_df[[1]]), ] # match the order of input

  # Rename column to matched the name of use_args
  if(!is.null(names(use_args))){
    # e.g., c(-x="label", "--yMin", "--yMax"), will rename the column 'label' in ref_df to 'x'
    use_args_vecName <- names(use_args)
    change_df <- use_args[use_args_vecName != ""] %>%
      tibble(use_arg = names(.), ref_col = .) %>%
      mutate(new_ref_colnames = str_remove(use_arg, "^\\-+"))
    # Change the ref_df column name to match
    ref_col_map <- dplyr::pull(change_df, new_ref_colnames, name = ref_col)
    colnames(ref_df)[colnames(ref_df) %in% names(ref_col_map)] <- ref_col_map[colnames(ref_df)[colnames(ref_df) %in% names(ref_col_map)]]

    # Change ref_arg names
    use_arg_map <- dplyr::pull(change_df, use_arg, name = ref_col)
    use_args[use_args_vecName != ""] <- use_arg_map[use_args[use_args_vecName != ""]]
  }
  # use_args_name <- str_remove(use_args, "^\\-+")

  # formatting arguments
  # Example use_args:
  # use_args <- c(-x="label", "--yMin", "--yMax")
  use_args_name <- str_remove(use_args, "^\\-+") # name of the column in ref_df
  arg_found_flag <- use_args_name %in% colnames(ref_df)
  if (any(!arg_found_flag)) {
    unfound_args <- use_args_name[!arg_found_flag]
    warning("the following arguments are unfound:\n\t", paste0(unfound_args, collapse = ", "))
  }
  use_args <- use_args[arg_found_flag]
  # use_args_name <- use_args_name[arg_found_flag]

  # Assembling arguments --------------------------------------------------------------------------
  args_list <- list()
  # arg_vec <- c()
  for (i in seq_along(use_args_name)) {
    cur_args_vals <- ref_df[[use_args_name[i]]]
    if ((length(unique(cur_args_vals)) == 1) & collapse_redundant) {
      cur_args_vals <- unique(cur_args_vals)
    }
    args_list[[use_args[i]]] <- cur_args_vals
    # arg_vec <- paste0(arg_vec, " ", use_args[i], " ", paste0(cur_args_vals, collapse = " "))
  }

  if(as_text){
    arg_vec <- ""
    for (i in seq_along(args_list)) {
      cur_arg <- args_list[[i]]
      if (length(cur_arg) == 0) {
        # Not sure if this will be the case
        arg_vec <- paste0(arg_vec, " ", use_args[i])
      } else {
        arg_vec <- paste0(arg_vec, " ", use_args[i], " ", paste0(cur_arg, collapse = " "))
      }
    }
    return(arg_vec)

  }else{
    names(args_list) <- str_remove(names(args_list), "^\\-+")
    return(args_list)
  }
}

# Command generator ===============================================================================
# MARK: ComputeMatrix Cmd
# A wrapper function to generate DeepTools' computeMatrix command =================================
DT_computeMatrix_cmd <- function(args_list, run_mode="reference-point", sample_table, sample_list=NULL,
                                 sample_path_colname="path", sample_label_colname="label", force = FALSE) {
  require(tibble)
  require(stringr)
  require(dplyr)
  # List known arguments for computeMatrix
  known_args = tribble(
    ~arg_name, ~flag_name, ~n_args, ~is_flag, ~mode, ~restrict_choice, ~required,
    "regions_file_name", c("--regionsFileName", "-R"), Inf, FALSE, "common", NA, TRUE,
    "score_file_name", c("--scoreFileName", "-S"), Inf, FALSE, "common", NA, TRUE,
    "out_file_name", c("--outFileName", "-o", "-out"), 1, FALSE, "common", NA, TRUE,
    "out_file_name_matrix", "--outFileNameMatrix", 1, FALSE, "common", NA, FALSE,
    "out_file_sorted_regions", "--outFileSortedRegions", 1, FALSE, "common", NA, FALSE,
    "reference_point", "--referencePoint", 1, FALSE, "reference-point", c("TSS", "TES", "center"), FALSE,
    "nan_after_end", "--nanAfterEnd", 0, TRUE, "reference-point", NA, FALSE,
    "region_body_length", c("--regionBodyLength", "-m"), 1, FALSE, "scale-regions", NA, FALSE,
    "start_label", "--startLabel", 1, FALSE, "scale-regions", NA, FALSE,
    "end_label", "--endLabel", 1, FALSE, "scale-regions", NA, FALSE,
    "unscaled_5prime", "--unscaled5prime", 0, TRUE, "scale-regions", NA, FALSE,
    "unscaled_3prime", "--unscaled3prime", 0, TRUE, "scale-regions", NA, FALSE,
    "before_region_start_length", c("--beforeRegionStartLength", "-b", "--upstream"), 1, FALSE, "common", NA, FALSE,
    "after_region_start_length", c("--afterRegionStartLength", "-a", "--downstream"), 1, FALSE, "common", NA, FALSE,
    "bin_size", c("--binSize", "-bs"), 1, FALSE, "common", NA, FALSE,
    "sort_regions", "--sortRegions", 1, FALSE, "common", c("descend", "ascend", "no", "keep"), FALSE,
    "sort_using", "--sortUsing", 1, FALSE, "common", c("mean", "median", "max", "min", "sum", "region_length"), FALSE,
    "sort_using_samples", "--sortUsingSamples", Inf, FALSE, "common", NA, FALSE,
    "average_type_bins", "--averageTypeBins", 1, FALSE, "common", c("mean", "median", "min", "max", "std", "sum"), FALSE,
    "missing_data_as_zero", "--missingDataAsZero", 0, TRUE, "common", NA, FALSE,
    "skip_zeros", "--skipZeros", 0, TRUE, "common", NA, FALSE,
    "min_threshold", "--minThreshold", 1, FALSE, "common", NA, FALSE,
    "max_threshold", "--maxThreshold", 1, FALSE, "common", NA, FALSE,
    "black_list_file_name", c("--blackListFileName", "-bl"), 1, FALSE, "common", NA, FALSE,
    "samples_label", "--samplesLabel", Inf, FALSE, "common", NA, FALSE,
    "smart_labels", "--smartLabels", 0, TRUE, "common", NA, FALSE,
    "quiet", "--quiet", 0, TRUE, "common", NA, FALSE,
    "scale", "--scale", 1, FALSE, "common", NA, FALSE,
    "verbose", "--verbose", 0, TRUE, "common", NA, FALSE,
    "number_of_processors", c("--numberOfProcessors", "-p"), 1, FALSE, "common", NA, FALSE,
    "metagene_center", "--metagene", 0, TRUE, "common", NA, FALSE,
    "transcript_id", "--transcriptID", 1, FALSE, "common", NA, FALSE,
    "exon_id", "--exonID", 1, FALSE, "common", NA, FALSE,
    "transcript_id_designator", "--transcript_id_designator", 1, FALSE, "common", NA, FALSE
  )
  if(!(run_mode %in% c("reference-point", "scale-regions"))){
    stop("Invalid run_mode")
  }
  known_args <- known_args[known_args$mode %in% c(run_mode, "common"), ]
  all_known_args <- unlist(known_args$flag_name)
  names(all_known_args) <- str_remove(all_known_args, "^\\-+")
  req_args <- known_args[known_args$required, ]
  req_args$flag_name <- lapply(req_args$flag_name, str_remove, pattern = "^\\-+")

  # Fetch sample args from sample table
  sample_args_list <- table_to_args(
    table_path = sample_table,
    use_entry = sample_list,
    # use_args = unlist(sample_args),
    use_args = c(`--samplesLabel` = sample_label_colname, `--scoreFileName` = sample_path_colname),
    collapse_redundant = TRUE
  )
  # if (F)  split_args_txt(sample_args, keep_dash = TRUE) # Internal check when debugging
  
  # Assemble args list
  A_list <- c(args_list, sample_args_list)
  A_list_fil <- A_list[names(A_list) %in% unique(names(all_known_args))]
  A_list_fil <- c(A_list_fil, sample_args_list)

  # Reorder A_list_fil to match order of all_known_args names
  reorder_idx <- match(names(all_known_args), names(A_list_fil))
  reorder_idx <- reorder_idx[!is.na(reorder_idx)]
  A_list_fil <- A_list_fil[reorder_idx]

  # Check for required arguments
  found_req_args <- sapply(req_args$flag_name, function(x) {
    any(x %in% names(A_list_fil))
  })
  if(!all(found_req_args)){
    missing_args <- req_args[!found_req_args, "arg_name"]
    stop(paste0("Missing required arguments: ", paste0(missing_args, collapse = ", ")))
  }

  # potential TODO: Check for redundant args?

  args_txt <- ""
  for (i in seq_along(A_list_fil)) {
    cur_arg <- A_list_fil[[i]]
    # Need to draw from all_known_args because we don't know the flag type of the current arg
    match_known_args <- all_known_args[names(all_known_args) %in% names(A_list_fil)[i]]
    if (is.logical(cur_arg)) {
      if (cur_arg) {
        args_txt <- paste0(args_txt, " ", match_known_args)
      }
    } else {
      args_txt <- paste0(args_txt, " ", match_known_args, " ", paste0(cur_arg, collapse = " "))
    }
  }
  
  cmd <- paste("computeMatrix", run_mode, args_txt)
  return(cmd)
}


# Main wrapper =====================================================================================
# MARK: Main wrapper

DT_wrapper <- function(run_what=c("computeMatrix"), args_list){
  cmd <- ""
  # computeMatrix ---------------------------------------------------------------------------------
  if ("computeMatrix" %in% run_what){
    # Fetch specific arguments
    # Set default
    run_mode <- "reference-point"
    use_samples <- NULL
    sample_path_col <- "path"
    sample_label_col <- "label"
    
    # Fetch args
    if ("run_mode" %in% names(args_list)) run_mode <- args_list$run_mode
    # Sample table - related args
    if ("sample_table" %in% names(args_list)) sample_table <- args_list$sample_table
    if("sample_path_col" %in% names(args_list)) sample_path_col <- args_list$sample_path_col
    if ("use_samples" %in% names(args_list)) use_samples <- args_list$use_samples
    
    cur_cmd <- DT_computeMatrix_cmd(args_list,
      run_mode = run_mode, sample_table = sample_table, sample_list = use_samples,
      sample_path_colname = "path", sample_label_colname = "label", force = FALSE
    )
    cmd <- paste0(c(cmd, cur_cmd), collapse = "\n")
  }
  return(cmd)
}




# Test section =====================================================================================
# MARK: Test section
if (FALSE) {
  # commandArgs()
  x1 <- c("/fs/pool/pool-r-library/packages/R-4.2.2-build_80xx/lib64/R/bin/exec/R", "--no-echo", "--no-restore", "--vanilla", "--file=/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/scripts/DT_wrapper.r", "--args", "-h", "--verbose", "--bw_info", "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/deeptools/bw_info.tsv", "--use_sample", "x", "y", "z", "asdf")
  # commandArgs(trailingOnly = TRUE)
  x2 <- c("-h", "--verbose", "--bw_info", "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/analysis/deeptools/bw_info.tsv", "--use_sample", "x", "y", "z", "asdf")
  x3 <- c("lead", "-h", "--verbose", "--bw_info", "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/analysis/deeptools/bw_info.tsv", "--use_sample", "x", "y", "z", "asdf")
  x4 <- c("lead", "lead2", "-h", "--verbose", "--bw_info", "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/analysis/deeptools/bw_info.tsv", "--use_sample", "x", "y", "z", "asdf", "--tail_flag")
  x5 <- c("DT_wrapper_2", "test", "--force", "-bed", "-bw_info", "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/analysis/deeptools/bw_info.tsv", "-bw", "Nr5a2_rep1", "Nr5a2_rep2", "Nr5a2_GSE218379_rep1", "Nr5a2_GSE218379_rep2", "Klf5_rep1", "Klf5_rep2", "ATAC_NT_RPKM", "ATAC_N_KD_RPKM", "ATAC_K_KD_RPKM", "ATAC_NK_DKD_RPKM", "H3K27ac_NT_RPKM", "H3K27ac_N_KD_RPKM", "H3K27ac_K_KD_RPKM", "H3K27ac_NK_DKD_RPKM", "--upstream", "1000", "--downstream", "1000", "--binSize", "50", "--kmeans", "5", "--clusterUsingSamples", "1", "6", "--refPointLabel", "center")
  x6 <- c(
    "DT_wrapper.r", "test",
    "--outFileName", "test_out.npz",
    "--bw_info", "/Volumes/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/analysis/deeptools/bw_info.tsv",
    "--referencePoint", "center", 
    "--regionsFileName", "test.bed",
    "--use_samples", "Nr5a2_rep1", "Nr5a2_rep2",
    "--asdf" # unrelated args
  )
  argv <- x6 # arg vector


  a_list <- format_argv(x2)
  a_list <- format_argv(x3)
  a_list <- format_argv(x4)
  a_list <- format_argv(x6, keep_dash = FALSE)

  a_list
  args_list <- a_list$arg
  use_args <- c("--colorMap", "--yMin")
  use_args <- c(`--sample`="label", "--path")
  args_list$bw_info = "/Volumes/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/analysis/deeptools/bw_info.tsv"
  table_to_args(table_path = args_list$bw_info, use_entry = NULL, use_args = use_args, collapse_redundant = TRUE)
  bw_info_path = args_list$bw_info
  # More real example
  a_list <- format_argv(x5)
  args_list <- a_list$arg
  arg_info <- a_list$info
  use_args <- c(`--bw` = "path", "--yMin", "--yMax", "--zMax", "--colorMap")
  table_to_args(
    table_path = args_list$bw_info,
    use_entry = args_list$bw,
    use_args = use_args,
    collapse_redundant = TRUE
  )
  DT_computeMatrix_cmd(args_list,
    run_mode = "reference-point",
    sample_table = args_list$bw_info, 
    sample_list = args_list$use_samples,
  )
  args_list
}



mode <- "test_retrieve_arg"
if (mode == "test_retrieve_arg") {
  # For getting test cmd args
  cat(commandArgs(), sep = '", "')

  cat("\n\n", rep("-", 20), "\n\n")
  cat(commandArgs(trailingOnly = TRUE), sep = '", "')
  cat("\n\n")

  cat("check bash param:\n")
  x <- system("( set -o posix ; set )", intern = T)
  print(subset(x, stringr::str_detect(x, "wd")))

  x2 <- system("env", intern = T)
  print(subset(x2, stringr::str_detect(x2, "wd")))

  x3 <- system("echo $wd", intern = T)
  print(x3)
}
