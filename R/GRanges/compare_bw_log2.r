#!/usr/bin/env Rscript 

compare_bw_log2 <- function(bw1, bw2, pseudocount=1, skip_zero_vs_zero=FALSE){
  bw_ovl <- findOverlaps(bw1, bw2)
  bw_ovl_df <- data.frame(b1_hit = queryHits(bw_ovl), b2_hit = subjectHits(bw_ovl))
  try({bw_ovl_df <- as_tibble(bw_ovl_df)})
  bw_ovl_df$b1_score <- bw1$score[bw_ovl_df$b1_hit]
  bw_ovl_df$b2_score <- bw2$score[bw_ovl_df$b2_hit]
  
  # Determine 0 vs. 0
  if(skip_zero_vs_zero){
    zero_v_zero_idx <- (bw_ovl_df$b1_score==0)&(bw_ovl_df$b2_score==0)
    bw_ovl_df_0 <- bw_ovl_df[zero_v_zero_idx , ]
    bw_ovl_df <- bw_ovl_df[!zero_v_zero_idx , ]
  }
  
  invisible(gc())
  # Calculate log2 FC
  bw_ovl_df$b1_score <- bw_ovl_df$b1_score + pseudocount
  bw_ovl_df$b2_score <- bw_ovl_df$b2_score + pseudocount
  
  bw_ovl_df$ratio <- bw_ovl_df$b1_score / bw_ovl_df$b2_score
  bw_ovl_df$log2FC <- log2(bw_ovl_df$ratio)
  
  # add the zero ones back
  if(skip_zero_vs_zero){
    bw_ovl_df_0$ratio <- 1
    bw_ovl_df_0$log2FC <- 0
    bw_ovl_df <- rbind(bw_ovl_df, bw_ovl_df_0)
    # Sort by b1, b2 order
    bw_ovl_df <- dplyr::arrange(bw_ovl_df, b1_hit, b2_hit)
    invisible(gc())
  }
  
  # Extract only the intersecting part
  # Assume that both bw1 and bw2 don't have self overlap regions (as they are bigwig format)
  bw_ovl_gr <- pintersect(bw1[bw_ovl_df$b1_hit], bw2[bw_ovl_df$b2_hit])
  invisible(gc())
  mcols(bw_ovl_gr) <- NULL
  bw_ovl_gr$score <- bw_ovl_df$log2FC
  bw_ovl_gr <- GenomicRanges::sort(bw_ovl_gr)
  
  return(bw_ovl_gr)
}


if(F){
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(GenomicRanges)
  library(rtracklayer)
  # library()
  # source_dir <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/dev/random_util/R"
  # source(paste0(source_dir, "/GRanges/annotate_overlap.r"))
  # source(paste0(source_dir, "/GRanges/findOverlaps_peaks.r"))
  # source(paste0(source_dir, "/get_loci_bw_score.r"))
  
  # Testing data set on ATAC seq of compound S treated cells
  data_dir <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/atac/merged/data"
  bw_dir <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/atac/merged/data/bw"
  outdir <- paste0(data_dir, "/bw_log2")
  bw1_path <- paste0(bw_dir, "/atac_2c_comS_mm10_wk91.127_bin1.spikeNorm_cellNorm.bw")
  bw2_path <- paste0(bw_dir, "/atac_2c_DMSO_mm10_wk90.128_bin1.spikeNorm_cellNorm.bw")
  
  
  bw1 <- BiocIO::import(bw1_path)
  bw2 <- BiocIO::import(bw2_path)
  gc()
  
  b1v2 <- compare_bw_log2(bw1, bw2, pseudocount=1, skip_zero_vs_zero=TRUE)
  gc()
  
  export.bw(b1v2, con=paste0(outdir, "/atac_2c_CmpS.wk91.127_vs_DMSO.wk90.128_bin1_log2FC.spikeNorm_cellNorm_manual.bw"))
  # # Binning genome
  # bs_mm10 <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  # bs_mm10_info <- as.data.frame(seqinfo(bs_mm10)) %>% 
  #   rownames_to_column("seqnames") %>% 
  #   dplyr::filter(str_detect(seqnames, "(^chr\\d+$)|(^chr[XYM]$)"))
  # 
  # bs_mm10_gr <- GRanges(seqnames=bs_mm10_info$seqnames, ranges=IRanges(start=1, end=bs_mm10_info$seqlengths))
  # bs_mm10_gr_bin10k <- unlist(tile(bs_mm10_gr, width=10E3))
  # bs_mm10_gr_bin1k <- unlist(tile(bs_mm10_gr, width=1E3))
}



# CLI =============================================================================================

if(!interactive() & sys.nframe()==0L){
  ## Not interactive and not being sourced from other function
  ## i.e., got run from command line interface.
  suppressPackageStartupMessages({
    library(argparser)
  })
  
  parser <- arg_parser("compare_bw_log2 | Calculate log2 fold-change of bw1 / bw2")
  
  parser <- add_argument(parser, short="-b1", arg="--bw1", 
                         help="Paths to input bigwig files, usually for test condition", 
                         type="character", default=NULL, nargs=1)
  
  parser <- add_argument(parser, short="-b2", arg="--bw2", 
                         help="Paths to input bigwig files, usually for control condition", 
                         type="character", default=NULL, nargs=1)
  
  parser <- add_argument(parser, arg="--pseudocount", help="Pesudocounts to add to both bw1 and bw2", type="numeric", default=1, nargs=1)
  
  # parser <- add_argument(parser, short="-p", arg="--n_core", 
  #                        help=paste0("Number of CPU for parallelization. Note the tool won't use more than one CPU for one input bigwig file."), 
  #                        type="numeric", default=1, nargs=1)
  
  parser <- add_argument(parser, short="-o", arg="--output", 
                         help=paste0("Path (preferably absolute) of the output bigwig file"), 
                         type="character", default=NULL, nargs=1)
  
  parser <- add_argument(parser, short="-z", arg="--skip_zero_vs_zero",
                         help="Skip comparison of region with zero value in both condition, this will automatically masked the log2 fold-change of these region as zero",
                         flag=TRUE)
  
  parser <- add_argument(parser, short="-O", arg="--overwrite",
                         help="Overwrite existing file",
                         flag=TRUE)
  
  parser <- add_argument(parser, short="-v", arg="--verbose",
                         help="verbosity",
                         flag=TRUE)
  
  try(argv <- parse_args(parser))
  verbose <- argv$verbose
  if(verbose){print(str(argv))}
  
  
  # Action time! ----------------------------------------------------------------------------------
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(GenomicRanges)
    library(rtracklayer)
    # library(parallel)
  })
  
  # Check if we should run this at all
  if(any(is.na(argv$output))){
    stop("Path to output wasn't provided")
  }
  
  outdir <- dirname(argv$output)
  if(!dir.exists(outdir)){
    stop("The output directory doesn't exist:\n\t", outdir, "\n")
  }
  
  if(tools::file_ext(argv$output) != "bw"){
    warning("The specificed output file doesn't contain .bw extension, adding to the output file name")
    argv$output <- paste0(argv$output, ".bw")
  }
  
  if(file.exists(argv$output)){
    if(!argv$overwrite){
      stop("Output file exist. Please use --overwrite if you would like to overwrite this file.")
    }else{
      if(verbose){cat("--overwrite specified; Removing exisiting output file.\n")}
      try({file.remove(argv$output)})
    }
  }
  
  # Read data
  if(verbose){cat("Reading input files ... ")}
  bw1 <- rtracklayer::import.bw(argv$bw1)
  bw2 <- rtracklayer::import.bw(argv$bw2)
  invisible(gc())
  if(verbose){cat("Done\n")}
  
  if(verbose){cat("Perform log2 fold-change of two bigwig files ... ")}
  bw_1v2 <- compare_bw_log2(bw1, bw2, pseudocount=argv$pseudocount, skip_zero_vs_zero=argv$skip_zero_vs_zero)
  if(verbose){cat("Done\n")}
  invisible(gc())
  # any(is.na(bw_1v2$score))
  
  if(verbose){cat("Saving output files ... ")}
  rtracklayer::export.bw(bw_1v2, con=argv$output)
  if(verbose){cat("Done\n")}
}

