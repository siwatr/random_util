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
