#!/usr/bin/env Rscript 

z_norm_bigwig <- function(gr, skip_zero=FALSE, NaN_as_zero=FALSE, ignore_seqnames=NULL, keep_original_score=FALSE){
  # Perform z-score normalization of bigwig data
  require(GenomicRanges)
  
  # Filter out unused bins
  # gr_fill will only be used for calculating standard deviation of the score
  gr_fil <- gr
  if(!is.null(ignore_seqnames)){
    gr_fil <- gr_fil[!tolower(seqnames(gr_fil)) %in% tolower(ignore_seqnames)]
  }
  
  # Dealing with unusual type of value
  nan_flag <- is.nan(gr_fil$score)
  na_flag <- is.na(gr_fil$score)
  inf_flag <- is.infinite(gr_fil$score)
  invalid_score_flag <- nan_flag | na_flag | inf_flag
  if(sum(invalid_score_flag) > 0){
    if(NaN_as_zero){
      gr_fil$score[invalid_score_flag] <- 0
    }else{
      # exclude from data set
      gr_fil <- gr_fil[!invalid_score_flag]
    }
  }
  
  if(skip_zero){
    gr_fil <- gr_fil[gr_fil$score!=0]
  }
  
  # calculate mean and standard diviation
  gr_fil_width <- width(gr_fil)
  mean_signal <- sum(gr_fil$score * gr_fil_width)/sum(gr_fil_width)
  sq_diff_signal <- (gr_fil$score - mean_signal)^2
  # var_signal <- mean((gr_fil$score - mean_signal)^2)
  var_signal <- sum(sq_diff_signal * gr_fil_width)/(sum(gr_fil_width)-1)
  sd_signal <- sqrt(var_signal)
  # hist(gr_fil$score, 100)
  
  # double check our work
  if(F){
    # For testing
    # gr_fil <- gr_fil[seqnames(gr_fil)=="chr1"]; gr_fil <- gr_fil[1:1000] # test
    # gr_fil <- gr_fil[seqnames(gr_fil)=="chr1"]; gr_fil <- gr_fil[(width(gr_fil)%%5)==0] # test
    
    # [DEPRECIATED] Calculate minimum bin
    # min_bin_length <- min(gr_fil_width, na.rm=TRUE)
    # good_bin_size <- all((gr_fil_width %% min_bin_length) == 0)
    # if(!good_bin_size){
    #   min_bin_length <- 1
    # }
    
    sig_vec <- sapply(seq_along(gr_fil), function(x){rep(gr_fil$score[x], each=gr_fil_width[x])})
    sig_vec <- c()
    for(i in seq_along(gr_fil)){
      sig_vec <- c(sig_vec, rep(gr_fil$score[i], each=gr_fil_width[i]))
    }
    length(sig_vec) == sum(gr_fil_width)
    mean(sig_vec); var(sig_vec); sd(sig_vec)
  }
  
  
  # Calculate z-score
  gr$z_score <- (gr$score - mean_signal)/sd_signal
  if(!keep_original_score){
    # Replace score field with z_score
    gr$score <- gr$z_score
    gr_mcols <- as.data.frame(mcols(gr))
    mcols(gr) <- gr_mcols[ , -which(colnames(gr_mcols)=="z_score"), drop=FALSE]
  }
  return(gr)
}

# Command line interface of the above function
z_norm_bigwig_cli <- function(bw_paths, outdir=".", out_suffix="_z_norm", 
                              skip_zero=FALSE, NaN_as_zero=FALSE, ignore_seqnames=NULL,
                              n_core=1, verbose=FALSE){
  require(rtracklayer)
  require(stringr)
  require(parallel)
  
  # In case user provide path with wildcard, attempt to extract one
  wildcard_check <- str_detect(bw_paths, "[*?]")
  if(any(wildcard_check)){
    if(verbose){cat("Interpreting wildcard character ... \n")}
    bw_paths_non_wc <- bw_paths[!wildcard_check]
    bw_path_wc <- bw_paths[wildcard_check]
    
    bw_path_wc_abs <- c() # absolute version of wildcard paths
    for(i in seq_along(bw_path_wc)){
      # bw_abs <- str_split_1(system(paste0("echo ", bw_path_wc[i]), intern=T), pattern="\\s+") # subject error when filename have space character in it.
      bw_abs <- system(paste0("readlink -f ", bw_path_wc[i]), intern=T)
      bw_path_wc_abs <- c(bw_path_wc_abs, bw_abs)
    }
    
    bw_paths <- c(bw_path_wc_abs, bw_paths_non_wc)
  }
  
  bw_check <- file.exists(bw_paths)
  if(all(!bw_check)){
    stop("None of the input files exist!")
  }else if(any(!bw_check)){
    missing_files <- paste0(bw_paths[!bw_check], collapse="\n")
    warning("The following files don't exist:\n", missing_files, "\n")
  }
  
  bw_paths <- bw_paths[bw_check]
  
  if(verbose){
    cat("\nFinal list of bw files ... \n")
    cat(bw_paths, sep="\n")
  }
  
  
  # Prepare paths for outputs
  bw_basenames <- tools::file_path_sans_ext(basename(bw_paths))
  out_paths <- paste0(outdir, "/", bw_basenames, out_suffix, ".bw")
  if(verbose){cat("\nZ-normalizing the bigwig files ... \n\n")}
  
  n_avail_core <- detectCores(logical=FALSE)
  n_core <- min(as.integer(abs(n_core)), n_avail_core, na.rm=TRUE)
  if(n_core == 1){
    bw_grl <- list()
    for(i in seq_along(bw_paths)){
      if(verbose){cat(bw_basenames[i], "\n")}
      bw <- rtracklayer::import.bw(bw_paths[i])
      bw_grl[[i]] <- z_norm_bigwig(bw, skip_zero=skip_zero, NaN_as_zero=NaN_as_zero, ignore_seqnames=ignore_seqnames)
    }
    
  }else{
    # Run normalization in parallel
    # NB: It seems that we can neither print something nor properly export bw files from inside mclapply function
    if(verbose){cat("Parallelly z-score normalize the bw files.\n")}
    out_list <- parallel::mclapply(seq_along(bw_paths), FUN=function(x){
      bw <- rtracklayer::import.bw(bw_paths[x])
      bw_z <- z_norm_bigwig(bw, skip_zero=skip_zero, NaN_as_zero=NaN_as_zero, ignore_seqnames=ignore_seqnames)
      
      tmp_list <- list(gr = bw_z, input = bw_paths[x])
      return(tmp_list)
    }, mc.cores=n_core)
    
    # Do this just because I'm paranoid that the order may get changed.
    input_vec <- sapply(out_list, function(x){x$input})
    # input_vec <- input_vec[c(4, 2, 1,3)]; match(input_vec, bw_paths) # Test: scramble the order
    bw_grl <- sapply(out_list, function(x){x$gr})
    x_order <- match(bw_paths, input_vec)
    input_vec <- input_vec[x_order]
    if(!all(input_vec == bw_paths)){
      warning("Something went wrong with matching the order of input. The saved files may not come from the correct input. Consider running the function with `n_core`=1")
    }
    
    bw_grl <- bw_grl[x_order]
  }
  
  if(F){
    # Testing parallization
    cl <- makeCluster(n_core, type="PSOCK")
    # initiation
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(GenomicRanges)
        library(rtracklayer)
        library(stringr)
      })
    }) %>% invisible()
    # Give some variable to child process
    clusterExport(cl, varlist=c("bw_paths", "bw_basenames", "out_paths"))
    # Doesn't seems to work
    clusterApplyLB(cl, x=seq_along(bw_paths), fun=function(x){
      print(bw_paths[x])
      bw <- rtracklayer::import.bw(bw_paths[x])
      rtracklayer::export.bw(bw, con=out_paths[x])
    })
    stopCluster(cl)
  }
  
  # Save files
  if(verbose){cat("Writing output files.\n")}
  for(i in seq_along(bw_grl)){
    if(verbose){cat("\t", basename(out_paths[i]), "\n", sep="")}
    rtracklayer::export.bw(object=bw_grl[[i]], con=out_paths[i])
  }
  
  # Do not return anything
}

if(F){
  bw_path <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/data/bw/Klf5_8C_rep1_WK155.filt.nodup.sort_bin1.RPKM.bw"
  bw <- import.bw(bw_path)
  gr <- bw
  
  # for multiple inputs and parallel testing
  bw_paths <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/data/bw/*RPKM*bw"
  skip_zero=T; NaN_as_zero=T; ignore_seqnames="chrM"; keep_original_score=FALSE
  outdir="/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/2_analysis/WK_Nr5a2_Klf5/data/bw_test"; out_suffix="_test"
  if(!dir.exists(outdir)){dir.create(outdir)}
  n_core=5
  verbose=T
}


# CLI =============================================================================================

if(!interactive() & sys.nframe()==0L){
  ## Not interactive and not being sourced from other function
  ## i.e., got run from command line interface.
  suppressPackageStartupMessages({
    library(argparser)
  })
  
  parser <- arg_parser("z_norm_bigwig.r | Z-score normalization of bigig files")
  
  parser <- add_argument(parser, short="-b", arg="--bw", 
                         help="Paths to input bigwig files to be z-score normalized", 
                         type="character", default=NULL, nargs=1)
  
  parser <- add_argument(parser, short="-p", arg="--n_core", 
                         help=paste0("Number of CPU for parallelization. Note the tool won't use more than one CPU for one input bigwig file."), 
                         type="numeric", default=1, nargs=1)
  
  parser <- add_argument(parser, short="-s", arg="--outfile_suffix", 
                         help=paste0("Suffix of the output file (without file extension)"), 
                         type="character", default="_z_norm", nargs=1)
  
  parser <- add_argument(parser, short="-o", arg="--outdir", 
                         help=paste0("Path (preferably absolute) to an output directory. "), 
                         type="character", default=".", nargs=1)
  
  parser <- add_argument(parser, short="-i", arg="--ignore_chr", 
                         help=paste0(
                           "character vector of chromosome name to be ignored for standard deviation calculation. ",
                           "Depending on your analysis, it can be useful to exclude some chromosome such as ChrM from sd calculation."
                         ), 
                         type="character", default="", nargs=Inf)
  
  
  parser <- add_argument(parser, short="-n", arg="--NaN_as_zero",
                         help="Treat NaN (Not-a-Number) values and other values of the kind (e.g., NA) as zero.",
                         flag=TRUE)
  
  parser <- add_argument(parser, short="-z", arg="--skip_zero",
                         help="exclude zero value from standard deviation calculation step",
                         flag=TRUE)
  
  parser <- add_argument(parser, short="-v", arg="--verbose",
                         help="verbosity",
                         flag=TRUE)
  
  try(argv <- parse_args(parser))
  if(argv$verbose){print(str(argv))}
  
  
  # Action time! ----------------------------------------------------------------------------------
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(parallel)
  })
  
  z_norm_bigwig_cli(bw_paths=argv$bw, 
                    outdir=argv$outdir, 
                    out_suffix=argv$outfile_suffix, 
                    skip_zero=argv$skip_zero, 
                    NaN_as_zero=argv$NaN_as_zero, 
                    ignore_seqnames=argv$ignore_chr,
                    n_core=argv$n_core, 
                    verbose=argv$verbose)
}
