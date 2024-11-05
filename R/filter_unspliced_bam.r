if(F) source("/fs/gpfs41/lv13/fileset03/pool/pool-toti-bioinfo/bioinfo/siwat_chad/lib/R/add_R_libPaths.r")

suppressPackageStartupMessages({
  library(dplyr, lib.loc="/fs/pool/pool-r-library/packages/R-4.2.2-build_80xx/lib64/R/library")
  library(tibble)
  library(tidyr)
  library(stringr)
  library(Rsamtools)
  library(GenomicAlignments)
  library(rtracklayer)
})
source("/fs/gpfs41/lv13/fileset03/pool/pool-toti-bioinfo/bioinfo/siwat_chad/dev/random_util/R/check_sam_flag.r")

test_bam_path <- "/fs/gpfs41/lv13/fileset03/pool/pool-toti-bioinfo/bioinfo/siwat_chad/dev/random_util/tmp/tmp_total_RNA_PE.bam"
gtf_path <- "/fs/gpfs41/lv13/fileset03/pool/pool-toti-bioinfo/bioinfo/siwat_chad/0_refData/genome/mus_musculus/mm10/annotated_genes/ensembl/GRCm38.102/Mus_musculus.GRCm38.102_withChr.gtf"

# Samtools-related params
samtools_path <- "/fs/pool/pool-totipotency-software/hpcl67/samtools/samtools-1.14/samtools"
samtools_sort_thread = 4
sort_output_by_position=TRUE

# params
bam_path <- test_bam_path
outfile_path <- paste0(tools::file_path_sans_ext(bam_path), "_unspliced_noExon.bam")
min_mapq <- 20
paired_end <- TRUE
# lib_strand <- "reverse"
lib_strand <- "forward"

# Exon overlap threshold
# exon_overlap_offset = 1
# Percent of the reads that should be part of the intron
min_percent_intron_overlap = 0
min_bp_intron_overlap = 5

verbose <- TRUE


bam_filter_nascent(bam_path=bam_path, 
                   outfile_path=paste0(tools::file_path_sans_ext(bam_path), "_unspliced_noExon.bam"), 
                   gtf_path=gtf_path, 
                   min_mapq=20, paired_end=TRUE, lib_strand="auto", 
                   min_percent_intron_overlap=0, min_bp_intron_overlap=5,
                   samtools_path=samtools_path, samtools_sort_thread=4, sort_output_by_position=TRUE, 
                   verbose=TRUE)


# Function section ================================================================================
bam_filter_nascent <- function(bam_path, outfile_path, gtf_path, min_mapq=20, paired_end=TRUE, lib_strand="auto", 
                               min_percent_intron_overlap=0, min_bp_intron_overlap=5,
                               samtools_path=NULL, samtools_sort_thread=4, sort_output_by_position=TRUE, verbose=FALSE){
  # Reading inputs
  # BAM
  if(verbose) cat("Reading BAM file ... ")
  bam_param <- ScanBamParam(what=c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "isize", "mrnm", "mate_status"))
  bam <- scanBam(bam_path, param=bam_param)
  bam <- as_tibble(do.call("DataFrame", bam))
  if(verbose) cat("Done\n")
  
  # bam_ga <- readGAlignmentPairs(file=bam_path, param=bam_param, use.names=TRUE)
  # Filtering by map quality
  if(verbose) cat("Read quality and gap (splicing) filtering ... ")
  if(!is.null(min_mapq)){
    bam <- dplyr::filter(bam, mapq >= min_mapq)
  }
  
  # paired-end specific
  if(paired_end){
    # mate mapping to the same ref sequence (rname)
    bam <- dplyr::filter(bam, rname == mrnm)
  }
  
  # ignore spliced read (containing N in the CIGAR string)
  bam <- dplyr::filter(bam, !str_detect(cigar, "\\d+N"))
  
  # select out any reads that its mate got filtered out
  bam <- dplyr::filter(bam, duplicated(qname) | duplicated(qname, fromLast=TRUE))
  
  # annotate pair orientation
  bam$is_R1 <- sam_flag_check(bam$flag, 64)
  # bam <- arrange(bam, qname, desc(is_R1))
  
  if(verbose) cat("Done\n")
  # Check full exon overlap -------------------------------------------------------------------------
  
  calculate_qwidth <- function(cigar_vec){
    # Test set
    if(F) cigar_vec <- c("110M", "3S79M4D28M", "46M1I43M20S", "2S108M", "110M")
    
    # Let's work with non-redundant set of CIGAR
    cigar_vec_2 <- unique(cigar_vec)
    # ignore the insertion and soft clipping
    cal_qwidth <- str_remove_all(cigar_vec_2, pattern="\\d+[IS]") %>% 
      str_extract_all(., pattern="\\d+(?=[:alpha:])") %>% 
      sapply(FUN=function(x){sum(as.numeric(x))})
    names(cal_qwidth) <- cigar_vec_2
    
    # Get length from the actual input
    o <- cal_qwidth[cigar_vec]
    if(any(is.na(o))){
      warning("NA is introduced in the output")
    }
    return(o)
  }
  
  # read GTF
  if(!is.null(gtf_path)){
    if(verbose) cat("Run exon filtering: \n")
    if(is.character(gtf_path)){
      gtf <- BiocIO::import(gtf_path)
      gtf <- gtf[gtf$type == "exon"]
      # ignore overlap exon position.
    }else{
      gtf <- gtf_path
    }
    gtf <- GenomicRanges::reduce(gtf)
    
    # Turning bam into GRange object
    # Calculating adjusted qwidth for each type of cigar string
    bam$adj_qwidth <- calculate_qwidth(bam$cigar)
    bam_gr <- GRanges(seqnames=bam$rname, 
                      ranges=IRanges(start=bam$pos, width=bam$adj_qwidth), 
                      strand=bam$strand,
                      qname = bam$qname, R1=bam$is_R1)
    
    # Automatically check library strandedness type (if specify)
    if(lib_strand=="auto"){
      rand_r1_idx <- sample(which(bam_gr$R1), min(length(bam_gr), 10E4))
      suppressWarnings({
        tmp_ovl <- GenomicRanges::findOverlaps(bam_gr[rand_r1_idx], gtf, ignore.strand=TRUE)
        strand_stats <- table(as.character(strand(gtf))[subjectHits(tmp_ovl)])
        strand_stats_pc <- strand_stats/sum(strand_stats)*100
        max_strand <- strand_stats_pc[which.max(strand_stats_pc)]
        if(max_strand>60){
          lib_strand <- c(`-`="reverse", `+`="forward")[names(max_strand)]
        }else{
          lib_strand <- "none"
        }
      })
      
      if(verbose){
        cat("Detected library strandedness: ", lib_strand, " (", signif(max_strand, 3), " of R1 reads)\n\n")
      }
    }
    
    # Check overlap based on library strand type
    rev_strand_map <- c(`-`="+", `+`="-", `.`="*", `*`="*")
    if(tolower(lib_strand) %in% c("reverse", "reversed", "rev")){
      # R2 already have the same orientation as the mapped genomic feature (e.g., genes)
      strand(bam_gr)[bam_gr$R1] <- rev_strand_map[as.character(strand(bam_gr)[bam_gr$R1])]
    }else if(tolower(lib_strand) %in% c("forward", "forwarded", "fwd")){
      # R1 already have the same orientation as the mapped genomic feature (e.g., genes)
      strand(bam_gr)[!bam_gr$R1] <- rev_strand_map[as.character(strand(bam_gr)[!bam_gr$R1])]
      
    }else{
      strand(bam_gr) <- "*"
    }
    
    # Find overlap with exon
    suppressWarnings({
      ovl <- findOverlaps(bam_gr, gtf, ignore.strand=FALSE)
      ovl_width <- width(pintersect(bam_gr[queryHits(ovl)], 
                                    gtf[subjectHits(ovl)], 
                                    ignore.strand=FALSE))
    })
    
    read_width <- width(bam_gr[queryHits(ovl)])
    
    ovl_df <- as_tibble(ovl)
    ovl_df$read_width <- width(bam_gr[queryHits(ovl)])
    ovl_df$exon_ovl_len <- ovl_width
    ovl_df$exon_strand <- as.character(strand(gtf))[subjectHits(ovl)]
    
    # Some reads may spand over multiple non redundant exons
    ovl_df$multi_ovl <- duplicated(ovl_df$queryHits) | duplicated(ovl_df$queryHits, fromLast=TRUE)
    if(any(ovl_df$multi_ovl)){
      dup_ovl <- ovl_df %>% 
        dplyr::filter(multi_ovl) %>% 
        group_by(queryHits, exon_strand) %>% 
        reframe(subjectHits=NA, 
                read_width=unique(read_width),
                exon_ovl_len=sum(exon_ovl_len),
                multi_ovl=TRUE)
      
      if(any(duplicated(dup_ovl$queryHits))){
        ## i.e., reads that overlap with exons from different strand
        # dup_ovl <- bind_rows(dup_ovl, mutate(dup_ovl, exon_ovl_len=exon_ovl_len-1)) # Test
        dup_ovl <- dup_ovl %>% 
          group_by(queryHits) %>% 
          mutate(exon_ovl_len=max(exon_ovl_len)) %>% 
          ungroup() %>% 
          unique()
      }
      
      ovl_df <- bind_rows(dplyr::filter(ovl_df, !multi_ovl), dup_ovl) %>% 
        arrange(queryHits)
    }
    
    # Calculate intron overlap length and filter out reads that overlap with exon too much
    excl_read_idx <- ovl_df %>%
      mutate(intron_ovl_len = read_width-exon_ovl_len,
             intron_ovl_pc = intron_ovl_len/read_width*100,
             intron_ovl_pass = (intron_ovl_len >= min_bp_intron_overlap) & (intron_ovl_pc >= min_percent_intron_overlap)) %>% 
      dplyr::filter(!intron_ovl_pass) %>% 
      dplyr::pull(queryHits) %>% 
      unique()
    
    bam$intron_ovl_pass <- TRUE
    bam$intron_ovl_pass[excl_read_idx] <- FALSE
    
    
    # keep a pair that at least one end isn't fully overlap with exon
    qname_2_keep <- dplyr::filter(bam, intron_ovl_pass) %>% 
      dplyr::pull(qname) %>% 
      unique()
    
    if(verbose){
      n_read_pair <- length(unique(bam$qname))
      cat("Keeping ", format(length(qname_2_keep), big.mark=","), 
          " out of ", format(n_read_pair, big.mark=","), " reads (", 
          signif(length(qname_2_keep)/n_read_pair*100, 3), "%)\n", sep="")
    }
    bam <- dplyr::filter(bam, qname %in% qname_2_keep)
  }
  
  if(paired_end){
    # Making sure each qname has duplicated (i.e., pair)
    bam <- dplyr::filter(bam, duplicated(qname) | duplicated(qname, fromLast=TRUE))
  }
  
  qname_2_keep <- unique(dplyr::pull(bam, qname))
  if(length(qname_2_keep)==0){
    stop("No reads left for filtering")
  }
  
  
  # using samtools to filter the reads --------------------------------------------------------------
  if(!is.null(samtools_path) & all(!is.na(samtools_path))){
    samtools <- samtools_path
  }else{
    # Assuming there is samtools available in user PATH variable
    samtools <- "samtools"
  }
  
  # export qname to be kept
  tmp_file_qname <- tempfile(tmpdir=dirname(bam_path), pattern=paste0("tmp_", basename(tools::file_path_sans_ext(bam_path)), "_qname_"), fileext=".txt")
  cat(qname_2_keep, file=tmp_file_qname, sep="\n", append=FALSE)
  
  # filter BAM file
  tmp_file_bam <- tempfile(tmpdir=dirname(outfile_path), pattern=paste0("tmp_unsorted_", basename(tools::file_path_sans_ext(outfile_path)), "_"), fileext=".bam")
  cmd <- paste0(samtools, " view -@ ", samtools_sort_thread, " --qname-file ", tmp_file_qname, " -o ", tmp_file_bam, " ", bam_path)
  if(verbose){cat(str_replace_all(cmd, "\\s+\\-", " \\\\ \n-"), "\n\n")}
  system(cmd)
  
  # sort bam
  if(verbose){cat("Sorting output ...\n")}
  if(sort_output_by_position){
    sort_param <- ""
  }else{
    sort_param <- " -n"
  }
  cmd <- paste0(samtools, " sort -@ ", samtools_sort_thread, sort_param, " -o ", outfile_path, " ", tmp_file_bam)
  if(verbose){cat(str_replace_all(cmd, "\\s+\\-", " \\\\ \n-"), "\n\n")}
  system(cmd)
  
  if(sort_output_by_position){
    cmd <- paste0(samtools, " index ", outfile_path)
    if(verbose){
      cat("Indexing BAM file ... \n")
      cat(str_replace_all(cmd, "\\s+\\-", " \\\\ \n-"), "\n\n")
    }
    system(cmd)
  }
  
  # Remove tmp files
  if(verbose){cat("clean up ...\n\n")}
  
  if(file.exists(outfile_path)){ # Check output file integrity before cleaning up
    # If there is some reads, remove temporary file
    n_reads_output <- system(paste0(samtools, " view -c ", outfile_path), intern=TRUE)
    n_reads_output <- as.double(n_reads_output)
    if(n_reads_output > 0){
      file.remove(tmp_file_qname)
      file.remove(tmp_file_bam)
    }else{
      warning("Cannot detect any reads in the output files, keep the temporary file for now.\n")
    }
  }
  
}



# tmp =============================================================================================
View(head(bam))
hist(bam$mapq, 100)
sum(bam$mapq<=20)

