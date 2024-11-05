# FIMO output header: 
# 1 motif_id
# 2 motif_alt_id
# 3 sequence_name
# 4 start
# 5 stop
# 6 strand
# 7 score
# 8 p-value
# 9 q-value
# 10 matched_sequence

fimo_2_bw <- function(fimo_path, outfile, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, relative_score=FALSE){
  fimo_df <- read.table(fimo_path, sep="\t", header=TRUE)
  fimo_gr <- makeGRangesFromDataFrame(dplyr::select(fimo_df, seqnames=sequence_name, start, stop, strand, score, p_value = p.value, q_value=q.value), keep.extra.columns=TRUE)
  # fimo_gr$score %>% hist(100)
  fimo_gr <- GenomicRanges::sort(fimo_gr)
  
  x <- findOverlaps(fimo_gr, fimo_gr, ignore.strand=T)
  fimo_self_ovl <- as_tibble(x) %>% 
    dplyr::filter(queryHits != subjectHits) %>% 
    nrow()
  
  if(fimo_self_ovl > 0){
    fimo_gr_red <- GenomicRanges::reduce(fimo_gr, ignore.strand=T)
    fimo_gr_red$score <- 0
    ovl <- findOverlaps(fimo_gr_red, fimo_gr, ignore.strand=T)
    # Inherit score from single overlap
    multi_ovl_q_idx <- queryHits(ovl)[duplicated(queryHits(ovl))]
    multi_ovl_flag <- queryHits(ovl) %in% multi_ovl_q_idx
    fimo_gr_red$score[queryHits(ovl)[!multi_ovl_flag]] <- fimo_gr$score[subjectHits(ovl)[!multi_ovl_flag]]
    
    multi_ovl_df <- as_tibble(ovl)[multi_ovl_flag , ]
    multi_ovl_df$s_score <- fimo_gr$score[multi_ovl_df$subjectHits]
    s_df <- multi_ovl_df %>% 
      group_by(queryHits) %>% 
      reframe(score = mean(s_score))
    
    fimo_gr_red$score[s_df$queryHits] <- s_df$score
    
    fimo_gr <- fimo_gr_red
  }
  
  if(relative_score){
    fimo_gr$score <- fimo_gr$score/max(fimo_gr$score, na.rm=T)
  }
  
  seqlengths(fimo_gr) <- seqlengths(genome)[names(seqlengths(fimo_gr))]
  
  # dir.exists(dirname(outfile))
  rtracklayer::export.bw(fimo_gr, con=outfile)
}




if(F){
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  # library()
  library(BiocIO)
  # "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/motifsearch/Nr5a2_mm10/fimo_MA0505.1_mm10/fimo.tsv"
  wd <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/motifsearch"
  
  # homeobox/OTX2: MA07122, NFYA: MA00603;TFAP2c MA08141
  fimo_path <- paste0(wd, "/Nr5a2_mm10/fimo_MA0505.1_mm10/fimo.tsv")
  outfile <- paste0(wd, "/Nr5a2_mm10/bw/Nr5a2_MA0505.1_mm10_score.bw")
  fimo_2_bw(fimo_path, outfile, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  
  fimo_path <- paste0(wd, "/JASPAR2022vertebrate_mm10/fimo/MA0599.1_KLF5/fimo.tsv")
  outfile <- paste0(wd, "/JASPAR2022vertebrate_mm10/bw/Klf5_MA0599.1_mm10_score.bw")
  fimo_2_bw(fimo_path, outfile, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  
  fimo_path <- paste0(wd, "/JASPAR2022vertebrate_mm10/fimo/MA0712.2_OTX2/fimo.tsv")
  outfile <- paste0(wd, "/JASPAR2022vertebrate_mm10/bw/Otx2_MA0712.2_mm10_score.bw")
  fimo_2_bw(fimo_path, outfile, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  
  fimo_path <- paste0(wd, "/JASPAR2022vertebrate_mm10/fimo/MA0060.3_NFYA/fimo.tsv")
  outfile <- paste0(wd, "/JASPAR2022vertebrate_mm10/bw/Nfya_MA0060.3_mm10_score.bw")
  fimo_2_bw(fimo_path, outfile, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  
  fimo_path <- paste0(wd, "/JASPAR2022vertebrate_mm10/fimo/MA0814.2_TFAP2C/fimo.tsv")
  outfile <- paste0(wd, "/JASPAR2022vertebrate_mm10/bw/TFAP2C_MA0814.2_mm10_score.bw")
  fimo_2_bw(fimo_path, outfile, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  
  fimo_path <- paste0(wd, "/jaspar_manual_download/mm10/Obox3_PH123.1/fimo/fimo.tsv")
  outfile <- paste0(wd, "/jaspar_manual_download/mm10/Obox3_PH123.1/bw/Obox3_PH123.1_mm10_score.bw")
  fimo_2_bw(fimo_path, outfile, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  
}
