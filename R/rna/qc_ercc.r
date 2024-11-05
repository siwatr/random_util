# Dependencies:
#  sce_extract.r

#' ERCC_QC_counts_fraction
#'
#' This function calculates the fraction of ERCC spike-in counts in a given count matrix.
#'
#' @param count_mat A matrix or a SummarizedExperiment object containing count data. 
#'                  If a SummarizedExperiment object is provided, the counts are extracted.
#' @param ercc_regex A regular expression pattern to identify ERCC spike-in rows. 
#'                   Default is "^ERCC".
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{sample}{Sample names (column names of the count matrix).}
#'   \item{spike_sum}{Sum of ERCC spike-in counts for each sample.}
#'   \item{total_sum}{Total sum of counts for each sample.}
#'   \item{pc_spike_in}{Percentage of ERCC spike-in counts relative to the total counts for each sample.}
#' }
#'
#' @details
#' The function first checks if the input is a SummarizedExperiment object and extracts the counts if necessary.
#' It then identifies the rows corresponding to ERCC spike-ins using the provided regular expression.
#' If no ERCC spike-ins are detected, the function stops with an error message.
#' Otherwise, it calculates the sum of ERCC spike-in counts and the total counts for each sample,
#' and returns a tibble with the calculated values and the percentage of spike-in counts.
#'
#' @examples
#' \dontrun{
#' # Example usage with a matrix
#' count_matrix <- matrix(data = rpois(100, lambda = 10), nrow = 10, ncol = 10)
#' rownames(count_matrix) <- c(paste0("ERCC_", 1:5), paste0("Gene_", 1:5))
#' result <- ERCC_QC_counts_fraction(count_matrix)
#'
#' # Example usage with a SummarizedExperiment object
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays = list(counts = count_matrix))
#' result <- ERCC_QC_counts_fraction(se)
#' }
#'
#' @importFrom SummarizedExperiment counts
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#' @export
ERCC_QC_counts_fraction <- function(count_mat, ercc_regex="^ERCC"){
  if(inherits(count_mat, "SummarizedExperiment")){
    count_mat <- counts(unsplitAltExps(count_mat, prefix.rows=F))
  }
  
  ercc_row_idx <- str_detect(rownames(count_mat), ercc_regex)
  if(!any(ercc_row_idx)){
    stop("No spike-in detected")
  }
  
  spike_in_df <- tibble(
    sample = colnames(count_mat),
    spike_sum = colSums(count_mat[ercc_row_idx , ]),
    total_sum = colSums(count_mat),
    pc_spike_in = (spike_sum/total_sum)*100
  )
  return(spike_in_df)
}


ERCC_extract <- function(sce, ercc_ref, assay_name="fpkm", tidy=FALSE){
  if(!inherits(sce, "SummarizedExperiment")){
    stop("The sce variable is expected to be a SummarizedExperiment object")
  }
  
  expr_df <- sce_extract(sce, query_regex="^ERCC", assay_names=assay_name, use_all_exp=TRUE)
  
  # Pick assays to show (the ERCC manual use FPKM)
  avail_assays <- unique(expr_df$assay)
  assay_name <- assay_name[assay_name %in% avail_assays]
  
  if(length(assay_name)>1){
    assay_name <- assay_name[1]
    warning("Expected the length of assay_name to be 1. Only the first assay ", assay_name, " will be used.")
  }
  
  if(length(assay_name)==0){
    assay_name <- avail_assays[1]
  }
  
  # expr_df is a tidy object by default
  expr_df <- dplyr::filter(expr_df, assay==assay_name)
  if(!tidy){
    expr_df <- tidyr::spread(expr_df, key="sample", value="value")
  }
  
  df <- full_join(expr_df, dplyr::rename(ercc_ref, feature_symbol=ercc_id), by="feature_symbol") %>% 
    tidyr::replace_na(list(value=0)) # TODO: add the data completion section after full_join
  
  # Rearrange the columns
  df <- dplyr::select(df, feature_symbol, any_of(colnames(ercc_ref)), any_of("assay"), everything())
  
  return(df)
}


ERCC_model <- function(ercc_df, sample_name=NULL, min_value=1, detected_only=FALSE){
  # min_counts: Minimum counts to consider each ERCC sequences to be detected.
  # detected_only: logical, if TRUE, perform linear model based on detected ERCC only.
  
  if(is.null(sample_name)){
    # cat(colnames(ercc_df), sep='", "')
    constant_ercc_df_colnames <- c("feature_symbol", "id", "subgroup", "conc_mix1", "conc_mix2", "expected_foldChange_ratio", "log2_mix1_mix2", "assay")
    sample_name = subset(colnames(ercc_df), !(colnames(ercc_df) %in% constant_ercc_df_colnames))
  }
  
  sample_found <- sample_name %in% colnames(ercc_df) # not the definitive way to filter for samples
  if(all(!sample_found)){
    stop("Couldn't find any samples")
  }else if(any(!sample_found)){
    n_missing <- sum(!sample_found)
    warning("Some samples are missing from ercc_df (", n_missing, "/", length(sample_found), "). ",
            "Missing samples:\n\t", paste0(sample_name[!sample_found], collapse=", "))
  }
  sample_name <- sample_name[sample_found]
  
  # How many ERCC sequence are detect?
  n_detected=colSums(ercc_df[ , sample_name] >= min_value)
  
  # Fitting model and check model's stats
  stats_df <- tibble()
  value_df <- tibble()
  for(i in seq_along(sample_name)){
    # Fitting linear model
    tmp_df <- ercc_df %>% 
      dplyr::select(feature_symbol, value=!!sym(sample_name[i]), conc_mix1, conc_mix2) %>% 
      mutate(detect = value >= min_value)
    
    tmp_df <- tmp_df %>% 
      mutate(original_idx = seq_along(feature_symbol)) %>% 
      arrange(desc(conc_mix1)) %>% 
      mutate(idx = seq_along(feature_symbol), 
             mx1_detect = idx <= max(idx[detect], na.rm=TRUE)) %>% 
      arrange(desc(conc_mix2)) %>% 
      mutate(idx = seq_along(feature_symbol), 
             mx2_detect = idx <= max(idx[detect], na.rm=TRUE)) %>% 
      arrange(original_idx) %>% 
      dplyr::select(-idx, -original_idx)
    
    model_df <- tmp_df %>% 
      mutate(log_mix1 = log10(conc_mix1), 
             log_mix2 = log10(conc_mix2), 
             log_value =log10(value + 0.001))
    
    if(detected_only){
      mx1_data = dplyr::filter(model_df, mx1_detect)
      mx2_data = dplyr::filter(model_df, mx2_detect)
    }else{
      mx1_data = mx2_data = model_df
    }
    
    mx1_lm <- lm(log_value ~ log_mix1, data=mx1_data)
    mx2_lm <- lm(log_value ~ log_mix2, data=mx2_data)
    # mx1_loess <- loess(log_value ~ log_mix1, data=model_df)
    # mx2_loess <- loess(log_value ~ log_mix2, data=model_df)
    
    # residuals(mx1_lm) %>% plot(seq_along(.), .)
    # coef(mx1_lm)
    model_df$lm_pred_log_mix1 <- predict(object=mx1_lm, newdata=model_df)
    model_df$lm_pred_log_mix2 <- predict(object=mx2_lm, newdata=model_df)
    # model_df$loess_pred_mix1 <- predict(object=mx1_loess, newdata=model_df)
    # model_df$loess_pred_mix2 <- predict(object=mx2_loess, newdata=model_df)
    
    value_df <- model_df %>% 
      mutate(sample=sample_name[i]) %>% 
      dplyr::select(feature_symbol, conc_mix1, conc_mix2, log_mix1, log_mix2, 
                    sample, value, log_value, 
                    lm_pred_log_mix1, lm_pred_log_mix2) %>% 
      bind_rows(value_df, .)
    
    stats_df <- tibble(
      sample=sample_name[i],
      lm_r2_mix1 = summary(mx1_lm)$r.squared,
      lm_r2_mix2 = summary(mx2_lm)$r.squared,
      lm_pcc_mix1 = cor(model_df$log_value, model_df$lm_pred_log_mix1, method="pearson"),
      lm_pcc_mix2 = cor(model_df$log_value, model_df$lm_pred_log_mix2, method="pearson"),
      lm_scc_mix1 = cor(model_df$log_value, model_df$lm_pred_log_mix1, method="spearman"),
      lm_scc_mix2 = cor(model_df$log_value, model_df$lm_pred_log_mix2, method="spearman"),
      # loess_mix1 = cor(model_df$log_value, model_df$loess_pred_mix1, method="spearman"),
      # loess_mix2 = cor(model_df$log_value, model_df$loess_pred_mix2, method="spearman")
      lm_coef_mix1 = coefficients(mx1_lm)[2],
      lm_coef_mix2 = coefficients(mx2_lm)[2]
    ) %>% 
      bind_rows(stats_df, .)
  }
  
  return(list(value = value_df, 
              stats = stats_df, 
              n_detected=n_detected))
}


ERCC_plot_scatter <- function(ercc_stats, sample_names=NULL, top=NULL, min_value=1, mix="all"){
  # top: integer, plot the top n worst and best fit samples
  df <- ercc_stats$value
  stats_df <- ercc_stats$stats
  ranked_samples <- stats_df %>% 
    arrange(desc(lm_r2_mix1)) %>% 
    dplyr::pull(sample)
  
  if(is.null(sample_names)){
    sample_names <- unique(df$sample)
  }
  
  if(!is.null(top)){
    top <- as.integer(top)
    n_samples <- length(ranked_samples)
    if(is.integer(top) & (top > 0)){
      sample_names <- c(subset(sample_names, sample_names %in% ranked_samples[1:top]),
                        subset(sample_names, sample_names %in% ranked_samples[(n_samples-top):n_samples])) %>% 
        unique()
      cat("using these samples: \n", paste0(sample_names, collapse=", "), sep="")
    }else{
      warning("Invalid `top` argument, only integer is accept")
    }
  }
  
  if(!all(mix %in% c(1,2))){
    mix = "all"
  }
  
  if(mix=="all"){
    mx_vec <- c(1, 2)
  }
  
  
  df <- dplyr::filter(df, sample %in% sample_names)
  stats_df <- dplyr::filter(stats_df, sample %in% sample_names)
  p_list <- list()
  for(i in seq_along(sample_names)){
    tmp_df <- dplyr::filter(df, sample == sample_names[i])
    s_df <- dplyr::filter(stats_df, sample == sample_names[i])
    
    # Determine which one is detected in each mix
    tmp_df <- tmp_df %>% 
      mutate(original_idx = seq_along(feature_symbol),
             detect=value >= min_value) %>% 
      arrange(desc(conc_mix1)) %>% 
      mutate(idx = seq_along(feature_symbol), 
             mx1_detect = idx <= max(idx[detect], na.rm=TRUE)) %>% 
      arrange(desc(conc_mix2)) %>% 
      mutate(idx = seq_along(feature_symbol), 
             mx2_detect = idx <= max(idx[detect], na.rm=TRUE)) %>% 
      arrange(original_idx) %>% 
      dplyr::select(-idx, -original_idx)
    
    p_mx <- list()
    for(m in mx_vec){
      # Extract current stats
      cur_R2 <- signif(s_df[ , paste0("lm_r2_mix", m), drop=TRUE], 2)
      cur_PCC <- signif(s_df[ , paste0("lm_pcc_mix", m), drop=TRUE], 2)
      cur_SCC <- signif(s_df[ , paste0("lm_scc_mix", m), drop=TRUE], 2)
      # cur_coef <- signif(s_df[ , paste0("lm_coef_mix", m), drop=TRUE], 2)
      subtitle_txt <- paste0("mix ", m, ": R2=", cur_R2, ", PCC=", cur_PCC, ", SCC=", cur_SCC)
      
      p_mx[[ m ]] <- tmp_df %>%
        mutate(label = if_else(!!sym(paste0("mx", m, "_detect")), "Within detection range", "Undetected")) %>% 
        ggplot(aes(x=!!sym(paste0("conc_mix", m))+0.01, y=value+0.01)) +
        geom_hline(yintercept=c(0+0.01), color="gray50") +
        geom_vline(xintercept=c(0+0.01), color="gray50") +
        geom_point(aes(shape=label), alpha=0.7, size=3) +
        scale_shape_manual(values=c(1, 19)) +
        scale_x_log10() + scale_y_log10() +
        geom_smooth(method="lm", formula=y~x, color="gray25") +
        labs(title=sample_names[i], subtitle=subtitle_txt,
             x="Known Conc. (mix1) + 0.01", y="Expression + 0.01") +
        theme_minimal() +
        # coord_fixed() +
        theme(legend.position="bottom")
    }
    
    if(length(p_mx) > 1){
      p <- cowplot::plot_grid(plotlist=p_mx, nrow=1)
    }else if(length(p_mx)==1){
      p <- p_mx[[1]]
    }else{
      p <- NULL
    }
    
    p_list[[ sample_names[i] ]] <- p
  }
  
  if(length(p_list)==1){
    p_list <- p_list[[1]]
  }
  
  return(p_list)
}

ERCC_plot_fit_quality <- function(ercc_stats){
  n_df <- data.frame(sample=names(ercc_stats$n_detected), n=ercc_stats$n_detected)
  p_n <- n_df %>%
    ggplot(aes(y=factor(sample, levels=rev(unique(sample))), x=n)) +
    geom_vline(xintercept=0, color="gray25") +
    geom_col(fill="gray25") +
    scale_x_continuous(limits=c(0,  92), breaks=c(seq(0, 80, by=10), 92)) +
    labs(x="Number of ERCC species detected", y="Samle Name") +
    theme_minimal()
  
  range_r2 <- range(c(ercc_stats$stats$lm_r2_mix1, ercc_stats$stats$lm_r2_mix2), na.rm=TRUE)
  p_r2 <- ercc_stats$stats %>% 
    ggplot(aes(x=lm_r2_mix1, y=lm_r2_mix2)) +
    geom_abline() +
    geom_point(alpha=0.7, size=3) +
    coord_fixed(xlim=range_r2, ylim=range_r2) +
    labs(title="R square (R2) of linear models between two ERCC mixes", 
         x="R2 of Mix1", y="R2 of Mix2") +
    theme_minimal()
  
  # guess which mix was used
  pc_mix1_better <- sum(ercc_stats$stats$lm_r2_mix1 > ercc_stats$stats$lm_r2_mix2)/nrow(ercc_stats$stats)*100
  guess_mix <- if_else(pc_mix1_better >= 0.5, 1, 2)
  r2_col <- paste0("lm_r2_mix", guess_mix)
  r2_vec <- ercc_stats$stats[ , r2_col, drop=TRUE]
  p_rank <- ercc_stats$stats %>%
    dplyr::select(sample, starts_with("lm_r2")) %>% 
    arrange(desc(!!sym(r2_col))) %>% 
    mutate(rank=seq_along(sample)) %>% 
    ggplot(aes(x=rank, y=!!sym(r2_col))) +
    geom_point(alpha=0.7, size=3) +
    labs(title="R square of each sample by rank", 
         x="Rranked Sample", y="R2") +
    ylim(c(min(0.5, r2_vec), 1)) +
    theme_minimal()
  
  out_list <- list(n_detect = p_n, 
                   R_square = p_r2,
                   R_square_rank = p_rank)
  
  return(out_list)
}

ERCC_QC <- function(sce, ercc_ref, assay_name="fpkm", outdir=NULL, outfile_prefix="", save_scatter_plots=TRUE){
  # A wrapper functions
  
  avail_assays <- assayNames(sce)
  if(!(assay_name %in% avail_assays)){
    warning("Couldn't find '", assay_name, "' in the input. Use '", avail_assays[1], "' instead. \nOther available assays are: \n\t", 
            paste0(avail_assays, collapse=", "))
    assay_name <- avail_assays[1]
  }
  
  ercc_df <- ERCC_extract(sce, ercc_ref, assay_name=assay_name, tidy=FALSE)
  ercc_stats <- ERCC_model(ercc_df, detected_only=F)
  ercc_scatter_plots <- ERCC_plot_scatter(ercc_stats)
  ercc_qc_plots <- ERCC_plot_fit_quality(ercc_stats)
  
  if(!is.null(outdir)){
    # Attempt to create the outdir
    if(!dir.exists(outdir)){
      dir.create(outdir, recursive=FALSE)
    }
    
    if(!dir.exists(outdir)){
      stop(outdir, " Doesn't exist!\n")
    }
    
    # Reformat prefix: make sure there is "_" at the end
    if(outfile_prefix != ""){
      outfile_prefix <- str_replace(paste0(outfile_prefix, "_"), "_+$", "_")
    }
    
    # Stats df
    n_df <- data.frame(sample = names(ercc_stats$n_detected), n_detect = ercc_stats$n_detected)
    left_join(ercc_stats$stats, n_df, by="sample") %>% 
      write.table(file=paste0(outdir, "/", outfile_prefix, "ERCC_QC_lm_stats.txt"), sep="\t", 
                  col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    if(save_scatter_plots){
      
    }
  }
}



ERCC_QC_old <- function(expr_df, ercc_ref, assay_name=character(0)){
  if(inherits(expr_df, "SummarizedExperiment")){
    expr_df <- sce_extract(expr_df, query_regex="^ERCC", assay_names="fpkm", use_all_exp=TRUE)
  }
  
  # Pick assays to show (the ERCC manual use FPKM)
  avail_assays <- unique(expr_df$assay)
  assay_name <- assay_name[assay_name %in% avail_assays]
  if(length(assay_name)==0){
    assay_name <- avail_assays[1]
  }
  
  expr_df <- dplyr::filter(expr_df, assay==assay_name)
  
  df <- full_join(expr_df, dplyr::rename(ercc_ref, feature_symbol=ercc_id), by="feature_symbol") %>% 
    tidyr::replace_na(list(value=0)) # TODO: add the data completion section after full_join
  
  tmpFn_model_cor <- function(df){
    # Fitting linear model
    model_df <- df %>% 
      dplyr::select(feature_symbol, sample, value, conc_mix1, conc_mix2) %>% 
      mutate(log_mix1 = log10(conc_mix1), log_mix2 = log10(conc_mix2), log_value =log10(value + 0.001))
    mx1_lm <- lm(log_value ~ log_mix1, data=model_df)
    mx2_lm <- lm(log_value ~ log_mix2, data=model_df)
    # mx1_loess <- loess(log_value ~ log_mix1, data=model_df)
    # mx2_loess <- loess(log_value ~ log_mix2, data=model_df)
    
    # residuals(mx1_lm) %>% plot(seq_along(.), .)
    # coef(mx1_lm)
    model_df$lm_pred_mix1 <- predict(object=mx1_lm, newdata=model_df)
    model_df$lm_pred_mix2 <- predict(object=mx2_lm, newdata=model_df)
    # model_df$loess_pred_mix1 <- predict(object=mx1_loess, newdata=model_df)
    # model_df$loess_pred_mix2 <- predict(object=mx2_loess, newdata=model_df)
    
    stats_df <- tibble(
      sample=unique(tmp_df$sample),
      lm_scc_mix1 = cor(model_df$log_value, model_df$lm_pred_mix1, method="spearman"),
      lm_scc_mix2 = cor(model_df$log_value, model_df$lm_pred_mix2, method="spearman"),
      # loess_mix1 = cor(model_df$log_value, model_df$loess_pred_mix1, method="spearman"),
      # loess_mix2 = cor(model_df$log_value, model_df$loess_pred_mix2, method="spearman")
      lm_coef_mix1 = coefficients(mx1_lm)[2],
      lm_coef_mix2 = coefficients(mx2_lm)[2]
    )
    
    return(list(value = model_df, cor = stats_df))
  }
  
  stats_df <- tibble()
  avail_sample <- unique(df$sample)
  p_list <- list()
  # avail_sample <- sort(sample(avail_sample, 5))
  for(s in avail_sample){
    tmp_df <- df %>% 
      dplyr::filter(sample==s) %>% 
      dplyr::arrange(conc_mix1) %>% 
      mutate(is_zero = value==0) 
    
    n_detected <- sum(!tmp_df$is_zero)
    pc_detected <- n_detected/nrow(tmp_df)*100
    
    # Concentration where we start to loose the detection power 
    min_undetected_idx <- min(which(!tmp_df$is_zero)) - 1
    tmp_df$detected <- TRUE
    if(min_undetected_idx > 0){
      tmp_df$detected[1:min_undetected_idx] <- FALSE
    }
    
    tmp_cor_df1 <- tmpFn_model_cor(tmp_df)$cor %>% 
      mutate(group="all")
    if(any(tmp_df$detected)){
      tmp_cor_df2 <- tmpFn_model_cor(df=dplyr::filter(tmp_df, detected))$cor %>% 
        mutate(group="detected")
    }else{
      tmp_cor_df2 <- tmp_cor_df1
      tmp_cor_df2[ , sapply(tmp_cor_df2, is.numeric)] <- 0
    }
    tmp_cor_df <- bind_rows(tmp_cor_df1, tmp_cor_df2) %>% 
      mutate(n_detect = n_detected,
             pc_detect = pc_detected)
    stats_df <- bind_rows(stats_df, tmp_cor_df)
    
    p_list[[s]] <- tmp_df %>% 
      gather(key="conc_set", "known_conc", conc_mix1, conc_mix2) %>% 
      mutate(value=value+0.001,
             conc_set = str_extract(conc_set, "mix\\d+")) %>% 
      ggplot(aes(x=known_conc, y=value)) +
      geom_smooth(formula=y~x, method="loess", color="gray25", alpha=0.1) +
      geom_smooth(formula=y~x, aes(group=detected), method="lm", color="orange", fill="darkgoldenrod1", alpha=0.1) +
      geom_line(alpha=0.5, color="gray80", linetype=2) +
      geom_vline(xintercept=0, linewidth=1, color="gray25") +
      geom_point(aes(group=detected, color=detected), alpha=0.3, size=4) +
      # geom_line() +
      facet_grid(.~conc_set) +
      scale_x_log10() + scale_y_log10() +
      labs(title=unique(tmp_df$sample), 
           subtitle=paste0("Detected: ", n_detected, "/", nrow(tmp_df), " (", signif(pc_detected, 3), " %)\n", 
                           "SCC (lm-all):              ", signif(tmp_cor_df1$lm_scc_mix1, 2), ", ", signif(tmp_cor_df1$lm_scc_mix2, 2), "\n",
                           "SCC (lm-detected):  ",        signif(tmp_cor_df2$lm_scc_mix1, 2), ", ", signif(tmp_cor_df2$lm_scc_mix2, 2)),
           x="Known Conc. (attomoles/ug)", 
           y=toupper(unique(tmp_df$assay))) +
      theme_minimal() +
      theme(legend.position="none", 
            axis.line.x=element_line(color="gray25", linewidth=1))
  }
  
  return(list(stats=stats_df,
              plots = p_list))
}



# Test ------------------------------------------------------------------------------------------

if(F){
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(SingleCellExperiment)
  library(DESeq2)
  
  library(ggplot2)
  library(ggbeeswarm)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(RColorBrewer)
  
  # Test input section ==============================================================================
  ercc_ref_path <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/0_refData/ERCC_spikein/ERCC_control_analysis_mod.txt"
  ercc_ref <- as_tibble(read.table(ercc_ref_path, sep="\t", header=TRUE, stringsAsFactors=FALSE))
  ercc_ref
  
  wd <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/rna_seq/AI/P829_20240712_Nr5a2_cKO_polyA"
  sce_path <- paste0(wd, "/results/star_featureCounts_cKO/star_cKO_featureCount_gene_body_gene_SCE.rds")
  # sce_path <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/rna_seq/AI/P776_20240321_Nr5a2Tev3/results/rsem_star_TevRef/RSEM_STAR_TevRef_filtered_gene_SCE.rds"
  sce_path <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/rna_seq/AI/P776_20240321_Nr5a2Tev3/results/20240716_matchExpression/rsem_star_TevRef/RSEM_STAR_TevRef_gene_SCE.rds"
  # sce_path <- paste0(wd, "/results/star_featureCounts/star_featureCounts_gene_body_gene_SCE.rds")
  sce_path <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/rna_seq/AI/P846_20240918_Nr5a2Tev/results/star_featureCounts_cKO/star_cKO_featureCount_gene_exon_gene_SCE.rds"
  sce <- readRDS(sce_path)
  
  # Test run ----------------------------------------------------------------------------------------
  assayNames(sce)
  sce_extract(sce, assay_names="counts", query_regex="^obox", search_field="gene_name", use_all_exp=TRUE, normalized=T, colData_colnames=c("condition", "cell_type", "bio_rep"), case_insensitive=T)
  expr_df <- sce_extract(sce, query_regex="^ERCC", assay_names="fpkm", use_all_exp=TRUE) # counts
  expr_df <- sce_extract(sce, assay_names="tpm", query_regex="^ERCC", use_all_exp=TRUE, colData_colnames=c("condition", "cell_type", "bio_rep"))
  expr_df <- sce_extract(sce, assay_names="tpm", query_regex="^ERCC", use_all_exp=TRUE, colData_colnames=c("condition", "cell_type", "bio_rep"), normalized=TRUE)
  ercc_counts <- ERCC_QC_counts_fraction(sce)
  qc_list <- ERCC_QC(expr_df, ercc_ref)
  ercc_m <- tpm(altExp(sce))
  ercc_stats_df <- full_join(qc_list$stats, ercc_counts, by="sample")
  
  # Actual test
  ercc_df <- ERCC_extract(sce, ercc_ref, assay_name="fpkm", tidy=FALSE)
  ercc_stats <- ERCC_model(ercc_df, detected_only=F)
  ercc_scatter_plots <- ERCC_plot_scatter(ercc_stats, top=5)
  ercc_scatter_plots
  ercc_qc_plots <- ERCC_plot_fit_quality(ercc_stats)
  
  
  ercc_stats$value
  ercc_stats$stats %>% 
    mutate(rank=rank(lm_r2_mix1)) %>% 
    arrange(rank) %>% 
    ggplot(aes(x=rank, y=lm_r2_mix1)) + 
    geom_point()
  
  ercc_stats$stats %>% 
    ggplot(aes(x=lm_r2_mix1, y=lm_r2_mix2)) +
    geom_abline() +
    geom_point() +
    theme_minimal() +
    coord_fixed()
  
  # correlation between two mix
  r_anno_vec <- dplyr::pull(ercc_ref, subgroup, name=ercc_id)
  r_anno_df <- ercc_ref %>% 
    mutate(log_conc_mix1 = log10(conc_mix1), log_conc_mix2 = log10(conc_mix2)) %>% 
    dplyr::select(ercc_id, subgroup, log_conc_mix1, log_conc_mix2) %>% 
    column_to_rownames(var="ercc_id")
  ref_col_fn <- colorRamp2(seq(-2, 5, length.out=8), viridis(8))
  ercc_m %>% 
    {.[names(r_anno_vec), ]} %>% 
    {log10(. + 0.001)} %>% 
    Heatmap(cluster_columns=F, cluster_rows=T, 
            col=viridis(8), 
            # right_annotation=rowAnnotation(sub_group = r_anno_vec),
            right_annotation=rowAnnotation(
              df=r_anno_df, 
              col=list(subgroup=setNames(cividis(length(unique(r_anno_df$subgroup))), 
                                         nm=unique(r_anno_df$subgroup)),
                       log_conc_mix1 = ref_col_fn, 
                       log_conc_mix2 = ref_col_fn)
            ),
            name="log_tpm")
  ercc_stats_df %>% 
    # dplyr::filter(group == "all") %>% 
    ggplot(aes(x=lm_scc_mix1, y=lm_scc_mix2)) +
    geom_abline() +
    geom_point() +
    coord_fixed() +
    facet_grid(.~group) +
    theme_minimal()
  
  # correlation between prediction correlation and linear model's coefficient
  ercc_stats_df %>% 
    dplyr::filter(group == "detected") %>%
    ggplot(aes(x=lm_scc_mix1, y=lm_coef_mix1)) +
    geom_smooth(method=lm) +
    # geom_hline(yintercept=0, linetype=2, color="orange") +
    geom_point() +
    # coord_fixed() +
    facet_grid(.~group) +
    theme_minimal()
  
  # Number of detected ERCC
  ercc_stats_df %>% 
    dplyr::select(sample, n_detect, pc_detect) %>% 
    unique() %>% 
    mutate(sample_factor = factor(sample, levels=rev(unique(sort(sample))))) %>% 
    ggplot(aes(x=n_detect, y=sample_factor)) +
    geom_col() +
    # xlim(c(0, 92)) +
    scale_x_continuous(breaks=c(seq(0, 90, by=10), 92), limits=c(0, 92)) +
    labs(x = "Number of detected ERCC species", y="Sample") +
    theme_minimal()
  
  # ERCC linear model prediction correlation
  ercc_stats_df %>% 
    gather(key="mix", value="scc", lm_scc_mix1, lm_scc_mix2) %>% 
    mutate(sample_factor = factor(sample, levels=rev(unique(sort(sample))))) %>% 
    ggplot(aes(x=scc, y=sample_factor, fill=group)) +
    # ggplot(aes(x=scc, y=sample, fill=group)) +
    geom_col() +
    geom_vline(xintercept=0, color="gray10") +
    facet_grid(group~mix) +
    theme_minimal() +
    labs(title="ERCC content QC", subtitle="Correlation between liner model predicted vs. observed ERCC spike-in expression", 
         x="Spearman's R", y="Sample") +
    scale_x_continuous(breaks=seq(-1,1, by=0.25), limits=c(NA,1)) +
    theme(legend.position="none")
  
  # Percent of spike-in reads
  ercc_stats_df %>% 
    dplyr::select(sample, pc_spike_in) %>% 
    unique() %>% 
    mutate(sample_factor = factor(sample, levels=rev(unique(sort(sample))))) %>% 
    ggplot(aes(x=pc_spike_in, y=sample_factor)) +
    geom_col() +
    labs(x="% Spike-in reads", y="Sample") +
    # xlim(c(0, max(ercc_stats_df$pc_spike_in))) +
    theme_minimal()
  
  # Number of ERCC vs. % ERCC in library
  ercc_stats_df %>% 
    dplyr::select(sample, n_detect, pc_spike_in) %>% 
    unique() %>% 
    ggplot(aes(x=n_detect, y=pc_spike_in)) +
    geom_hline(yintercept=0, color="gray10") +
    geom_smooth(method="lm", color="orange", fill="darkgoldenrod1", alpha=0.2) +
    geom_point(size=3, alpha=0.4) +
    labs(x="Number of detected ERCC species", y="% ERCC read counts") +
    theme_minimal()
  
  col_data <- as.data.frame(colData(sce)) %>% 
    dplyr::select(-any_of("sample")) %>% 
    rownames_to_column(var="sample")
  
  require(ggbeeswarm)
  ercc_stats_df %>% 
    dplyr::select(-any_of("group")) %>% unique() %>% 
    left_join(., col_data, by="sample") %>% 
    ggplot(aes(x = condition, y=pc_spike_in, color=condition)) +
    geom_boxplot(width=0.5) +
    geom_quasirandom(size=5, alpha=0.7) +
    # geom_beeswarm(size=5, cex=5, alpha=0.7) +
    ylim(c(0, NA)) +
    labs(title="Percent of spike-in content", x = "Condition", y="% Spike-in reads") +
    theme_minimal()
  
  ercc_stats_df %>% 
    ggplot(aes(x=lm_scc_mix1)) +
    geom_histogram(bins=100) +
    xlim(c(0,1)) +
    facet_grid(.~group) +
    theme_minimal()
  
  # worst sample
  s1_1 <- ercc_stats_df %>% 
    dplyr::filter(group=="all") %>% 
    arrange(lm_scc_mix1) %>% 
    dplyr::pull(sample) %>% {.[1]}
  s1_2 <- ercc_stats_df %>% 
    dplyr::filter(group=="detected") %>% 
    arrange(lm_scc_mix1) %>% 
    dplyr::pull(sample) %>% {.[1]}
  # best sample
  s2_1 <- ercc_stats_df %>% 
    dplyr::filter(group=="all") %>% 
    arrange(desc(lm_scc_mix1)) %>% 
    dplyr::pull(sample) %>% {.[1]}
  s2_2 <- ercc_stats_df %>% 
    dplyr::filter(group=="detected") %>% 
    arrange(desc(lm_scc_mix1)) %>% 
    dplyr::pull(sample) %>% {.[1]}
  
  
  p_list[[s1_1]]
  p_list[[s1_2]]
  p_list[[s2_1]]
  p_list[[s2_2]]
  
}


# CLI =============================================================================================

