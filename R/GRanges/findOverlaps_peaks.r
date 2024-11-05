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

