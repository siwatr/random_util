get_loci_bw_score <- function(query, bw, query_refPoint="exact", query_extend=0, value_type="all", out_col_prefix="bw_"){
  require(dplyr)
  require(IRanges)
  require(GenomicRanges)
  require(rtracklayer)
  # Read an average genome coverage score
  # query = GRanges object of a BED file
  # bw = path to a bigwig file or GRanges object of a bigwig file
  # query_refPoint, query_extend: query modification parameters:
  #   query_extend: <numeric>: modifying the query region: default = NA, meaning the exact input region will be used in this function.
  #   query_refPoint %in% c("exact", "center"): whether the extension will be done from the center or edge of the input region
  # method %in% c("mean", "median): mode data representation, default: mean.
  # value_type %in% c("all", "sum", "mean", "sd", "median", "max")
  
  ## Checking input files: ------------------------------------------------------------------------
  # query <- check_GR(query)
  if(inherits(bw, "GRanges")){
    bw_score <- bw
  }else if(is.character(bw) & file.exists(bw)){
    bw_score <- rtracklayer::import.bw(bw)
  }else{
    stop("incorrect bw input: incorrect format")
  }
  
  if(!("score" %in% colnames(mcols(bw_score)))){
    stop("incorrect bw input: no score field")
  }
  
  ## Modify the query regions ---------------------------------------------------------------------
  # Default setting: query_refPoint="exact"; query_extend=NA
  # Test setting: query_refPoint="center"; query_extend=100
  # Verify setting
  if(query_extend>0){
    query_refPoint <- tolower(query_refPoint)
    if(!(query_refPoint %in% c("exact", "center"))){
      warning("Incorrect `query_refPoint` argument, skip query modification step!")
    }else{
      query <- switch(query_refPoint,
                      exact = query,
                      center = resize(query, width=1, fix="center"))
      
      query <- flank(query, width=query_extend, both=TRUE)
    }
  }
  
  ## Check value type input variable 
  accept_val <- c("all", "sum", "mean", "sd", "median", "max")
  if(!all(value_type %in% accept_val)){
    warning(paste0("`value_type` contain unaccepted value. The function can only accept the following value: ", paste0(accept_val, collapse=", ")))
    
    if(any(accept_val %in% value_type)){
      value_type <- accept_val[accept_val %in% value_type]
    }else{
      # Revert to default value
      value_type <- "all"
    }
  }
  
  if("all" %in% value_type){
    ## Use all values, except "all"
    use_value_type <- accept_val[-1]
  }else{
    use_value_type <- value_type
  }
  
  ## Acquiring coverage data: ---------------------------------------------------------------------
  hit <- findOverlaps(query, bw_score)
  hit_df <- data.frame(qHit = queryHits(hit),
                       sHit = subjectHits(hit)) %>%
    mutate(qWidth = width(query[qHit]),
           in_width = width(pintersect(query[qHit], bw_score[sHit])), # intersect width: The real overlap regions between the two
           score = bw_score[sHit]$score, ## Score of that bin
           sum_score = score * in_width) ## Area under the bar plot of the intersect region
  
  ## NB: the intersect area (in_width) doesn't include the bw region where the score = 0.
  ##      This is needed to be accounted for when calculating sd and median values.
  ## For each qHit (i.e., each query region), find score.
  hit_sum <- hit_df %>%
    group_by(qHit) %>%
    reframe(qWidth = mean(qWidth), # Should give the same number anyway
            zero_score_width=qWidth - sum(in_width), # qWidth outside the intersect region
            ## Calculate scores
            
            score_max = if("max" %in% use_value_type){max(score)}else{NA}, # i.e., the peak height
            score_sum = sum(sum_score),                                    # Area under the curve. This is NOT the read count
            score_mean = if("mean" %in% use_value_type){score_sum/qWidth}else{NA},
            ## calculation of SD and median might take longer time because of the rep() function (?)
            ## The calculation of SD doesn't accounted for the case where all the value is uniform, which might cause the problem
            score_sd = try(if("sd" %in% use_value_type){sd(c(numeric(zero_score_width), rep(score, in_width)))}else{NA}),
            score_median = if("median" %in% use_value_type){median(c(numeric(zero_score_width), rep(score, in_width)))}else{NA}) 
  
  ## add these scores to bw 
  query_mcols <- as_tibble(mcols(query))
  ## Add a new empty column to mcols
  score_colname <- c(paste0(out_col_prefix, "score_", use_value_type))
  query_mcols[ , score_colname] <- 0
  query_mcols[hit_sum$qHit , score_colname] <- hit_sum[ , paste0("score_", use_value_type)]
  mcols(query) <- query_mcols
  
  # ## For manually check the result when debugging the code:
  # tmp_df <- filter(hit_df, qHit==1875)
  # zero_vals <- numeric(mean(tmp_df$qWidth) - sum(tmp_df$in_width))
  # score_vec <- rep(tmp_df$score, tmp_df$in_width)
  # sum(c(zero_vals, score_vec))
  # median(c(zero_vals, score_vec))
  # sd(c(zero_vals, score_vec))
  # mean(c(zero_vals, score_vec))
  # max(zero_vals, score_vec)
  
  return(query)
}


# # Test version ====================================================================================
# 
# get_loci_bw_score <- function(query, bw, query_refPoint="exact", query_extend=0, value_type="all", out_col_prefix="bw_"){
#   require(dplyr)
#   require(IRanges)
#   require(GenomicRanges)
#   require(rtracklayer)
#   # Read an average genome coverage score
#   # query = GRanges object of a BED file
#   # bw = path to a bigwig file or GRanges object of a bigwig file
#   # query_refPoint, query_extend: query modification parameters:
#   #   query_extend: <numeric>: modifying the query region: default = NA, meaning the exact input region will be used in this function.
#   #   query_refPoint %in% c("exact", "center"): whether the extension will be done from the center or edge of the input region
#   # method %in% c("mean", "median): mode data representation, default: mean.
#   # value_type %in% c("all", "sum", "mean", "sd", "median", "max")
#   
#   ## Checking input files: ------------------------------------------------------------------------
#   # query <- check_GR(query)
#   if(inherits(bw, "GRanges")){
#     bw_score <- bw
#   }else if(is.character(bw) & file.exists(bw)){
#     bw_score <- rtracklayer::import.bw(bw)
#   }else{
#     stop("incorrect bw input: incorrect format")
#   }
#   
#   if(!("score" %in% colnames(mcols(bw_score)))){
#     stop("incorrect bw input: no score field")
#   }
#   
#   ## Modify the query regions ---------------------------------------------------------------------
#   # Default setting: query_refPoint="exact"; query_extend=NA
#   # Test setting: query_refPoint="center"; query_extend=100
#   # Verify setting
#   if(query_extend>0){
#     query_refPoint <- tolower(query_refPoint)
#     if(!(query_refPoint %in% c("exact", "center"))){
#       warning("Incorrect `query_refPoint` argument, skip query modification step!")
#     }else{
#       query <- switch(query_refPoint,
#                       exact = query,
#                       center = resize(query, width=1, fix="center"))
#       
#       query <- flank(query, width=query_extend, both=TRUE)
#     }
#   }
#   
#   ## Check value type input variable 
#   accept_val <- c("all", "sum", "mean", "sd", "median", "max")
#   if(!all(value_type %in% accept_val)){
#     warning(paste0("`value_type` contain unaccepted value. The function can only accept the following value: ", paste0(accept_val, collapse=", ")))
#     
#     if(any(accept_val %in% value_type)){
#       value_type <- accept_val[accept_val %in% value_type]
#     }else{
#       # Revert to default value
#       value_type <- "all"
#     }
#   }
#   
#   if("all" %in% value_type){
#     ## Use all values, except "all"
#     use_value_type <- accept_val[-1]
#   }else{
#     use_value_type <- value_type
#   }
#   
#   ## Acquiring coverage data: ---------------------------------------------------------------------
#   hit <- findOverlaps(query, bw_score)
#   hit_df <- data.frame(qHit = queryHits(hit),
#                        sHit = subjectHits(hit)) %>%
#     mutate(qWidth = width(query[qHit]),
#            in_width = width(pintersect(query[qHit], bw_score[sHit])), # intersect width: The real overlap regions between the two
#            score = bw_score[sHit]$score, ## Score of that bin
#            sum_score = score * in_width) ## Area under the bar plot of the intersect region
#   
#   ## NB: the intersect area (in_width) doesn't include the bw region where the score = 0.
#   ##      This is needed to be accounted for when calculating sd and median values.
#   ## For each qHit (i.e., each query region), find score.
#   hit_sum <- hit_df %>%
#     group_by(qHit) %>%
#     reframe(qWidth = mean(qWidth), # Should give the same number anyway
#             zero_score_width=qWidth - sum(in_width), # qWidth outside the intersect region
#             ## Calculate scores
#             
#             score_max = if("max" %in% use_value_type){max(score)}else{NA}, # i.e., the peak height
#             score_sum = sum(sum_score),                                    # Area under the curve. This is NOT the read count
#             score_mean = if("mean" %in% use_value_type){score_sum/qWidth}else{NA},
#             ## calculation of SD and median might take longer time because of the rep() function (?)
#             ## The calculation of SD doesn't accounted for the case where all the value is uniform, which might cause the problem
#             score_sd = try(if("sd" %in% use_value_type){sd(c(numeric(zero_score_width), rep(score, in_width)))}else{NA}),
#             score_median = if("median" %in% use_value_type){median(c(numeric(zero_score_width), rep(score, in_width)))}else{NA}) 
#   
#   ## add these scores to bw 
#   query_mcols <- as_tibble(mcols(query))
#   ## Add a new empty column to mcols
#   score_colname <- c(paste0(out_col_prefix, "score_", use_value_type))
#   query_mcols[ , score_colname] <- 0
#   query_mcols[hit_sum$qHit , score_colname] <- hit_sum[ , paste0("score_", use_value_type)]
#   mcols(query) <- query_mcols
#   
#   # ## For manually check the result when debugging the code:
#   # tmp_df <- filter(hit_df, qHit==1875)
#   # zero_vals <- numeric(mean(tmp_df$qWidth) - sum(tmp_df$in_width))
#   # score_vec <- rep(tmp_df$score, tmp_df$in_width)
#   # sum(c(zero_vals, score_vec))
#   # median(c(zero_vals, score_vec))
#   # sd(c(zero_vals, score_vec))
#   # mean(c(zero_vals, score_vec))
#   # max(zero_vals, score_vec)
#   
#   return(query)
# }
