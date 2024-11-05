#!/usr/bin/env Rscript 

# Main location: /fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/dev/random_util/R/GRanges

# Utility functions ===============================================================================
set_GRangesList_names <- function(peak_obj, peak_names=NULL, prefix="p_", show_warning=FALSE){
  # peak_obj: GRangesList or simple vector (e.g., to peak paths)
  warn_msg <- ""
  if(!is.null(peak_names)){
    # i.e., peak_names is provided: but is it correct?
    if(any(duplicated(peak_names))){
      warn_msg <- "Duplicated peak_names detected. "
      peak_names=NULL
    }
    if(length(peak_names) != length(peak_obj)){
      warn_msg <- "The lenght of peak_names is not the same as input peak_obj "
      peak_names=NULL
    }
  }
  
  # Attempt 1: use names associated with peak_obj
  if(is.null(peak_names) & !is.null(names(peak_obj))){
    p_names <- names(peak_obj)
    if(!any(duplicated(p_names))){
      peak_names <- p_names
    }else{
      warn_msg <- paste0(warn_msg, "Name of the peak_obj has some duplication. ")
    }
  }
  
  # Attempt 2: use file name
  if(is.null(peak_names) & inherits(peak_obj, "character")){
    # assume paths is provided
    p_base <- tools::file_path_sans_ext(basename(peak_obj))
    warn_msg <- paste(warn_msg, "Attempt to use file basename as peak name. ")
    if(!any(duplicated(p_base))){
      peak_names <- p_base
    }else{
      warn_msg <- paste0(warn_msg, "File basename is duplicated. Using the entry index as peak name instead. ")
    }
  }
  
  # Last Attempt x: use peak index
  if(is.null(peak_names)){
    warn_msg <- paste0(warn_msg, "File basename is duplicated. Using the entry index as peak name instead. ")
    peak_names <- paste0(prefix, str_pad(seq_along(peak_obj), width=str_length(length(peak_obj)), pad="0", side="left"))
  }
  if(show_warning){warning(warn_msg)}
  
  return(peak_names)
}

filter_peak_width <- function(peak_list, min_width=50, max_width=1E3, verbose=FALSE){
  GRanges_input <- inherits(peak_list, "GRanges")
  if(GRanges_input){
    peak_list <- GRangesList(peak_list)
    names(peak_list) <- "GRanges_1"
  }
  
  if(inherits(peak_list, "list")){
    peak_list <- GRangesList(peak_list)
    try(names(peak_list) <- set_GRangesList_names(peak_list, peak_names=names(peak_list), prefix="peak_"))
  }
  
  if(!inherits(peak_list, "CompressedGRangesList")){
    stop("Incorrect input format: ", class(peak_list))
  }
  
  # Filter peaks
  if(verbose){cat("Filtering peaks by width ...\n")}
  for(i in seq_along(peak_list)){
    w_vec <- width(peak_list[[i]])
    size_pass_flag <- (w_vec>=min_width) & (w_vec<=max_width)
    
    if(verbose){
      pc_keep <- signif(sum(size_pass_flag)/length(peak_list[[i]])*100, 3)
      cat("\t", names(peak_list)[i], ":\t", 
          format(sum(size_pass_flag), nsmall=2, big.mark=","), "/", 
          format(length(peak_list[[i]]), nsmall=2, big.mark=","), 
          " (", pc_keep, " %)\n", sep="")
    }
    # filter
    peak_list[[i]] <- peak_list[[i]][size_pass_flag]
  }
  
  if(GRanges_input){
    peak_list <- unlist(peak_list)
  }
  
  return(peak_list)
}


read_peak_list <- function(peak_paths, peak_names=NULL, min_width=0, max_width=Inf, verbose=FALSE){
  names(peak_paths) <- set_GRangesList_names(peak_paths, peak_names=peak_names, prefix="peak_")
  
  peak_list <- GRangesList()
  for(i in seq_along(peak_paths)){
    peak_list[[ names(peak_paths)[i] ]]  <- BiocIO::import(peak_paths[i])
  }
  
  # Check if we need to filter the peak by size
  size_mat <- sapply(peak_list, function(x){range(width(x))})
  do_size_filter <- any(size_mat[1, , drop=TRUE] < min_width) | any(size_mat[2, , drop=TRUE] > max_width)
  if(do_size_filter){
    peak_list <- filter_peak_width(peak_list, min_width=min_width, max_width=max_width, verbose=verbose)
  }
  
  return(peak_list)
}

# Peak overlap ====================================================================================

findOverlaps_peaks <- function(gr1, gr2, min_overlap_bp=1, min_overlap_percent=0, ignore.strand=FALSE){
  ovl <- findOverlaps(gr1, gr2, ignore.strand=ignore.strand)
  q_hit <- queryHits(ovl)
  s_hit <- subjectHits(ovl)
  
  ovl_width <- width(pintersect(gr1[q_hit], gr2[s_hit], ignore.strand=ignore.strand))
  ovl_pc <- (ovl_width / pmin(width(gr1[q_hit]), width(gr2[s_hit])))*100
  
  # filter overlaps
  ovl_filter_1 <- ovl_width >= min_overlap_bp
  ovl_filter_2 <- ovl_pc >= min_overlap_percent
  ovl_fil <- ovl[ovl_filter_1 & ovl_filter_2]
  
  return(ovl_fil)
}

annotate_overlap <- function(query_gr_list, subject_gr_list=NULL, ovl_col_prefix="",
                             min_overlap_bp=1, min_overlap_percent=0, ignore.strand=FALSE,
                             reduce_query=FALSE){
  # require(glue)
  if(is.null(subject_gr_list)){
    # i.e., self overlap check
    subject_gr_list <- query_gr_list
  }
  
  if(reduce_query){
    query_gr_list <- GRangesList(GenomicRanges::reduce(unlist(query_gr_list)))
    # print(length(query_gr_list)) # debug
  }
  
  for(q in seq_along(query_gr_list)){
    q_gr <- query_gr_list[[q]]
    # q_name <- names(query_gr_list)[q]
    
    for(s in seq_along(subject_gr_list)){
      s_gr <- subject_gr_list[[s]]
      s_name <- names(subject_gr_list)[s]
      ovl_col <- paste0(ovl_col_prefix, s_name)
      
      # Check overlap
      ovl <- findOverlaps_peaks(gr1=q_gr, gr2=s_gr, 
                                min_overlap_bp=min_overlap_bp, 
                                min_overlap_percent=min_overlap_percent, 
                                ignore.strand=ignore.strand)
      q_hit <- queryHits(ovl)
      df <- as.data.frame(mcols(q_gr))
      df[[ ovl_col ]] <- FALSE
      df[[ ovl_col ]][unique(q_hit)] <- TRUE
      
      # Update the current query granges
      mcols(q_gr) <- df
    }
    
    # Update the grange list
    mcols(query_gr_list[[q]]) <- mcols(q_gr)
  }
  
  if(reduce_query){
    query_gr_list <- unlist(query_gr_list)
  }
  
  return(query_gr_list)
}


macs2_reproducible_peaks <- function(peak_list, 
                                     min_rep_ovl=1, # minimum number of replicate with that peak
                                     min_width=50, max_width=1E3,
                                     min_overlap_bp=50, min_overlap_percent=50, 
                                     ignore.strand=TRUE, verbose=FALSE){
  
  # Peak size filtering
  if(is.null(names(peak_list))){
    names(peak_list) <- paste0("p_", str_pad(seq_along(peak_list), width=str_length(length(peak_list)), pad="0", side="left"))
  }
  
  size_mat <- sapply(peak_list, function(x){range(width(x))})
  do_size_filter <- any(size_mat[1, , drop=TRUE] < min_width) | any(size_mat[2, , drop=TRUE] > max_width)
  if(do_size_filter){
    peak_list <- filter_peak_width(peak_list, min_width=min_width, max_width=max_width, verbose=verbose)
  }
  
  # Looking for overlap
  anno_peak <- annotate_overlap(peak_list, 
                                ovl_col_prefix="peak_ovl_", 
                                min_overlap_bp=min_overlap_bp, 
                                min_overlap_percent=min_overlap_percent, 
                                reduce_query=TRUE, # will collapse the entired GRangesList into a single GRange object
                                ignore.strand=ignore.strand)
  
  if(inherits(anno_peak, "list") | inherits(anno_peak, "CompressedGRangesList")){
    if(length(anno_peak)>1){warning("Expecting anno_peak to be a list of one element.")}
    anno_peak <- unlist(anno_peak)
  }
  
  # width(anno_peak) %>% hist(100) # debug
  df <- as.data.frame(mcols(anno_peak)) %>% 
    dplyr::select(starts_with("peak_ovl"))
  anno_peak$n_rep_ovl <- rowSums(df)
  
  # Filter overlap
  anno_peak <- anno_peak[anno_peak$n_rep_ovl >= min_rep_ovl]
  return(anno_peak)
}


macs2_reproducible_peaks_cli <- function(peak_paths, peak_names=NULL,
                                         min_rep_ovl=NULL, # minimum number of replicate with that peak
                                         min_width=50, max_width=1E3,
                                         min_overlap_bp=50, min_overlap_percent=50, 
                                         ignore.strand=TRUE,
                                         outfile_path=NULL, overwrite=FALSE,
                                         do_return=FALSE,
                                         verbose=FALSE){
  # A wrapper to work with file path inputs
  peak_list <- read_peak_list(peak_paths=peak_paths, peak_names=peak_names, 
                              min_width=min_width, max_width=max_width, verbose=verbose)
  if(is.null(min_rep_ovl)){
    # i.e., need to be in all reps
    min_rep_ovl <- length(peak_list)
  }
  
  # print(peak_list)
  # verbose option is for peak filtering
  rep_peak <- macs2_reproducible_peaks(peak_list=peak_list, 
                                       min_width=min_width, max_width=max_width, 
                                       min_rep_ovl=min_rep_ovl, min_overlap_bp=min_overlap_bp, 
                                       min_overlap_percent=min_overlap_percent, verbose=verbose)
  
  # Determine whether we should write the output, or should we print to stdout?
  do_write_outfile <- FALSE
  if(!is.null(outfile_path)){
    if(!file.exists(outfile_path) | overwrite){
      do_write_outfile <- TRUE
    }else if(file.exists(outfile_path)){
      warning("Output path is already exisit, print the output instead")
    }
  }
  
  # Writing output to file
  if(do_write_outfile){
    ext <- tools::file_ext(outfile_path)
    known_ext=c("bed", "rds", "xlsx", "csv", "tsv", "txt")
    if(!(ext %in% known_ext)){
      ext <- "bed" # set default to bed file for IGV viewing
    }
    
    if(tolower(ext)=="bed"){
      export.bed(rep_peak, con=outfile_path)
      
    }else if(tolower(ext)=="rds"){
      saveRDS(rep_peak, outfile_path)
      
    }else if(tolower(ext)=="xlsx"){
      writexl::write_xlsx(as.data.frame(rep_peak), path=outfile_path)
      
    }else if(tolower(ext)=="csv"){
      write.csv(as.data.frame(rep_peak), file=outfile_path, row.names=FALSE, col.names=TRUE, quote=FALSE)
      
    }else if(tolower(ext) %in% "tsv" | "txt"){
      write.table(as.data.frame(rep_peak), file=outfile_path, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
    
  }else if(do_return){
    # useful for debugging
    return(rep_peak)
    
  }else{
    # i.e., print out to stdout: This might not be the best option for R
    df <- as.data.frame(rep_peak)
    df <- df[ , c("seqnames", "start", "end", "strand")]
    write.table(df, sep="\t", file=stdout(), row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
}


# Testing -----------------------------------------------------------------------------------------
if(F){
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(stringr)
  
  script_dir="/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/dev/random_util/R/GRanges"
  source(paste0(script_dir, "/findOverlaps_peaks.r"))
  source(paste0(script_dir, "/annotate_overlap.r"))
  
  peak_dir="/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/data/peak/narrow"
  peak_basenames=c(
    N_2C_rep1="Nr5a2_2C_untreated_rep1_WK178_peaks.narrowPeak",
    N_2C_rep2="Nr5a2_2C_untreated_rep2_WK181_peaks.narrowPeak",
    N_4C_rep1="Nr5a2_4C_untreated_rep2_WK182_peaks.narrowPeak"
  )
  peak_paths <- file.path(peak_dir, peak_basenames)
  names(peak_paths) <- names(peak_basenames)
  peak_list <- GRangesList()
  for(i in seq_along(peak_paths)){
    p_name <- names(peak_paths)[i]
    gr <- BiocIO::import(peak_paths[i])
    gr <- gr[(width(gr)>=50) & (width(gr)<=1E3)]
    peak_list[[p_name]] <- gr
  }
  
  peak_list
  
  anno_peak_list <- annotate_overlap(peak_list, ovl_col_prefix="ovl_", min_overlap_bp=50, min_overlap_percent=50, ignore.strand=TRUE, reduce_query=T)
  width(anno_peak_list) %>% hist(100)
  df <- as.data.frame(mcols(anno_peak_list))
  
  
  rep_peak <- macs2_reproducible_peaks_cli(peak_paths=peak_paths, do_return=TRUE)
}


# CLI =============================================================================================
if(!interactive() & sys.nframe()==0L){
  ## Not interactive and not being sourced from other function
  ## i.e., got run from command line interface.
  suppressPackageStartupMessages({
    library(argparser)
  })
  
  parser <- arg_parser("macs2_reproducible_peaks.r | Determin the reproducible peaks from multiple input files")
  
  parser <- add_argument(parser, short="-i", arg="--input", 
                         help="Paths to input files in BED format, or outputs from MACS2 peak calling (e.g., narrowPeak, or broadPeak files)", 
                         type="character", default=NULL, nargs=Inf)
  
  parser <- add_argument(parser, short="-o", arg="--out", 
                         help=paste0("Path (preferably absolute) to an output files. Currently support the following formats: ",
                                     "bed, rds, xlsx, csv, tsv, and txt. If not provided, the ouptut will be printed out the stdout."), 
                         type="character", default=NULL, nargs=1)
  
  parser <- add_argument(parser, short="-n", arg="--peak_names",
                         help=paste0("Names of the peaks associated with each peak files. ",
                                     "The peak names have to be unique, and have the same number of elements as --input option. ",
                                     "Note that this won't affect the BED or Stdout output"),
                         type="character", default=NULL, nargs=Inf)
  
  parser <- add_argument(parser, short="-R", arg="--min_rep",
                         help=paste0("Minimum number of replicates (each peak files) that have to be part of the overlap to consider it overlap. ",
                                     "By default, all replicates need to be part of the overlap of each reproducible peaks"),
                         type="numeric", default=NULL, nargs=1)
  
  parser <- add_argument(parser, short="-m", arg="--min_width",
                         help=paste0("Minimum width (in bp) of peaks"),
                         type="numeric", default=50, nargs=1)
  
  parser <- add_argument(parser, short="-X", arg="--max_width",
                         help=paste0("Maximum width (in bp) of peaks"),
                         type="numeric", default=1000, nargs=1)
  
  # Overlap params
  parser <- add_argument(parser, short="-B", arg="--min_overlap_bp",
                         help=paste0("Minimum width (in bp) of overlap"),
                         type="numeric", default=50, nargs=1)
  
  parser <- add_argument(parser, short="-P", arg="--min_overlap_percent",
                         help=paste0("Minimum percentage of a smaller peaks that take part into the overlap. ",
                                     "For example, if we consider overlap peaks of width 50, and 100 bp, ",
                                     "if their overlap is 5 bp, then the overlap percentage of this pair will be 5/50*100 = 10 % of the smaller peaks"),
                         type="numeric", default=50, nargs=1)
  
  parser <- add_argument(parser, short="-S", arg="--stranded",
                         help="The overlap need to taking strand of the peaks/regions into account",
                         flag=TRUE)
  
  parser <- add_argument(parser, short="-O", arg="--overwrite",
                         help="Overwrite existing files",
                         flag=TRUE)
  
  parser <- add_argument(parser, short="-v", arg="--verbose",
                         help="verbosity",
                         flag=TRUE)
  
  try(argv <- parse_args(parser))
  if(argv$verbose){print(str(argv))}
  
  
  # RUN -------------------------------------------------------------------------------------------
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    # library(parallel)
    library(dplyr)
    library(stringr)
    library(tools)
  })
  
  # format some arguments
  if(any(is.na(argv$input))){
    stop("No input file provided!")
  }
  
  # Peak names
  if(any(is.na(argv$peak_names))){
    peak_names = NULL # would this work?
  }else if(length(argv$peak_names) != length(argv$input)){
    peak_names = NULL
  }else{
    peak_names = argv$peak_names
  }
  
  # min rep ovl
  if(any(is.na(argv$min_rep))){
    min_rep_ovl = NULL
  }else{
    min_rep_ovl = argv$min_rep
  }
  
  # Output path
  if(any(is.na(argv$out))){
    outfile_path = NULL
  }else{
    outfile_path = argv$out
  }
  
  # RUN
  macs2_reproducible_peaks_cli(
    peak_paths = argv$input, 
    # additional peak reading params
    peak_names = peak_names,
    min_rep_ovl = min_rep_ovl, # minimum number of replicate with that peak
    min_width=argv$min_width, 
    max_width=argv$max_width,
    min_overlap_bp=argv$min_overlap_bp, 
    min_overlap_percent=argv$min_overlap_percent, 
    ignore.strand=!argv$stranded,
    # output params
    outfile_path=outfile_path, 
    overwrite=argv$overwrite,
    do_return=FALSE
  )
}



