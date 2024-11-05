extract_SAMflag <- function(f){
  # Reference individual flag bit according to the SAM file specification
  # we'll check from bigger bit to the smaller one, hence the decreasing sorting
  flag_bits <- sort(2^c(0:11), decreasing=TRUE)
  
  # in case user provide the exact sam flag bit
  if(f %in% flag_bits){return(f)}
  
  # Otherwise it has to be break down
  f_vec <- c() # All bit information
  for(i in seq_along(flag_bits)){
    cur_diff <- f - flag_bits[i]
    if(cur_diff >= 0){
      # yes it has this flag bit
      f_vec <- c(f_vec, flag_bits[i])
      f <- cur_diff
    }
  }
  
  # for sanity check: f must be completely broken down to zero
  if(f != 0){stop("Impossible SAM flag is given (?)")}
  
  return(f_vec)
}



sam_flag_check <- function(query_flag, ref_flag, debug_mode=FALSE){
  # Check if query_flag (e.g., flag of individual read in sam/bam file) contain the information of the ref_flag (e.g., flag used for -f or -F in samtools view command)
  # e.g., query flag 163 contains flag bit 128
  
  check_single_flag <- function(q, r, r_info){
    if(q == r){
      return(TRUE)
    }else if(q < r){
      return(FALSE)
    }else{
      # r_info <- extract_SAMflag(r)
      q_info <- extract_SAMflag(q)
      return(all(r_info %in% q_info))
    }
  }
  
  # This version would be more likely to work when query_flag is given as an interger vector
  
  # Making a dictionery of flag inclusion
  r_info <- extract_SAMflag(ref_flag)
  q_unq <- sort(unique(query_flag))
  flag_map <- c()
  for(i in seq_along(q_unq)){
    flag_map[as.character(q_unq)[i]] <- check_single_flag(q=q_unq[i], r=ref_flag, r_info=r_info)
  }
  check <- flag_map[as.character(query_flag)]
  
  ## older version that use sapply
  # r_info <- extract_SAMflag(ref_flag)
  # check <- sapply(query_flag, FUN=function(q){check_single_flag(q=q, r=ref_flag, r_info=r_info)})
  
  if(debug_mode){
    message("flag mapping:\n")
    print(flag_map)
    cat("\n\n")
  }
  
  return(check)
}

# Test data ---------------------------------------------------------------------------------------
if(F){
  f <- query_flag <- 163
  query_flag <- c(83, 163, 147, 83)
  ref_flag <- 2
  ref_flag <- 64
  
  extract_SAMflag(163)
  sam_flag_check(query_flag=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ref_flag=c(2))
}

