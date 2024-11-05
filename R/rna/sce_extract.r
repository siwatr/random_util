# Calculate quality control matrix of dataset with ERCC spike-in data
# https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_086340.pdf

# Function section ================================================================================
# sce_2_expr_tidy: --------------------------------------------------------------------------------
#' Extract expression data (assay) and sample data (colData) from SingleCellExperiment object, 
#' then convert it to Data Frame in Tidy format
#'
#' This function converts a SingleCellExperiment object into a tidy data frame
#' with optional normalization and selection of specific assays and column data.
#'
#' @param sce A SingleCellExperiment object containing the single-cell RNA-seq data.
#' @param assay_names A character vector specifying the names of assays to include in the output. 
#'                    If empty, all assays will be included.
#' @param colData_colnames A character vector specifying the column names from colData to include in the output.
#'                         If empty, all columns will be included. If FALSE, no colData will be included.
#' @param normalized A logical value indicating whether to normalize the expression values using size factors.
#'                   Default is FALSE.
#'
#' @return A tibble containing the tidy data frame with columns for feature symbols, sample names, assay names,
#'         expression values, and selected colData columns.
#'
#' @details
#' The function performs the following steps:
#' 1. Converts the colData of the SingleCellExperiment object to a data frame and selects specified columns.
#' 2. Filters and ensures the inclusion of the "sample" column in the colData.
#' 3. Selects specified assays from the SingleCellExperiment object.
#' 4. Normalizes the expression values if the `normalized` parameter is TRUE and size factors are available.
#' 5. Extracts expression values from the assays and converts them to a tidy format.
#' 6. Merges the expression data with the selected colData and returns the result.
#'
#' @examples
#' # Example usage:
#' # sce <- SingleCellExperiment::SingleCellExperiment(...)
#' # tidy_df <- sce_2_expr_tidy(sce, assay_names = c("counts", "logcounts"), colData_colnames = c("cell_type", "condition"), normalized = TRUE)
#'
#' @importFrom dplyr select any_of left_join mutate bind_rows
#' @importFrom tibble as_tibble rownames_to_column tibble
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment colData assays assay
#' @importFrom SingleCellExperiment sizeFactors
#' @importFrom stats setNames
#' @export
sce_2_expr_tidy <- function(sce, assay_names=character(0), colData_colnames=character(0), normalized=FALSE){
  df <- as.data.frame(colData(sce)) %>% 
    dplyr::select(-any_of("sample")) %>% 
    rownames_to_column(var="sample") %>% 
    as_tibble()
  
  # Filter colData columns
  if(is.logical(colData_colnames)){
    if(!colData_colnames){
      # i.e., no colData will be included in the output
      colData_colnames <- "sample"
    }
  }
  
  colData_colnames <- colData_colnames[colData_colnames %in% colnames(df)] # pre-filter
  colData_colnames <- unique(c("sample", colData_colnames)) # makesure we have sample in there
  if(length(colData_colnames)==0){
    colData_colnames <- colnames(df)
  }
  df <- df[ , colData_colnames]
  
  # Make tidy
  assay_names <- assay_names[assay_names %in% names(assays(sce))]
  if(length(assay_names)==0){
    assay_names <- names(assays(sce))
  }
  
  if(is.null(sizeFactors(sce)) & normalized){
    warning("Couldn't detect pre-calculated size factor of the input, skip normalization.\n")
    normalized <- FALSE
  }
  
  if(normalized){
    sf <- tibble(sample=colnames(sce), factor=sizeFactors(sce))
  }
  
  # Extract experssion value
  expr_df <- tibble()
  for(i in seq_along(assay_names)){
    d <- as.data.frame(assay(sce, assay_names[i])) %>% 
      rownames_to_column(var="feature_symbol") %>% 
      as_tibble() %>% 
      gather(key="sample", value="value", -feature_symbol) %>% 
      mutate(assay=assay_names[i]) %>% 
      dplyr::select(feature_symbol, sample, assay, value, everything())
    
    if(normalized){
      d <- left_join(d, sf, by="sample") %>% 
        mutate(value=value/factor) %>% 
        dplyr::select(-factor)
    }
    
    expr_df <- bind_rows(expr_df, d)
  }
  
  # Merge to colData and return
  if(ncol(df)>1){
    out_df <- left_join(expr_df, df, by="sample", suffix=c("", ".y"))
  }else{
    out_df <- expr_df
  }
  
  return(out_df)
}


# sce_extract
#' Extract and Tidy Expression Data for Specific Genes from a SummarizedExperiment Object
#'
#' This function extracts expression data for genes matching a specified regular expression 
#' from a SummarizedExperiment object and returns it in a tidy format.
#'
#' @param sce A SummarizedExperiment object. This can be an object of class SummarizedExperiment, 
#' DESeqDataSet, or SingleCellExperiment.
#' @param query_regex A character string specifying the regular expression to match gene names. 
#' Default is "^ERCC".
#' @param search_field A character string specifying which columns of rowData should be used for searching
#' Default is "rownames", which will use the rownames of the SummarizedExperiment object as the search reference.
#' @param assay_names A character string specifying the assay to extract data from. Default is "counts".
#' @param colData_colnames A character vector specifying the column names from colData to include in the output. 
#' Default is an empty character vector.
#' @param use_all_exp A logical value indicating whether to use all experiments in the SummarizedExperiment object. 
#' Default is FALSE.
#' @param exp_name A character string specifying the name of the alternative experiment to use. 
#' If NULL, the main experiment is used. Default is NULL.
#' @param normalized A logical value indicating whether to return normalized expression values. Default is FALSE.
#'
#' @return A tidy data frame containing the expression data for the specified genes.
#'
#' @examples
#' # Assuming `sce` is a SummarizedExperiment object
#' expr_df <- sce_extract(sce, query_regex="^ERCC", assay_names="counts", colData_colnames=c("sample_id"), 
#'                        use_all_exp=FALSE, exp_name=NULL, normalized=FALSE)
#'
#' @import SummarizedExperiment
#' @import stringr
#' @export

sce_query_gene <- function(sce, query_regex="^ERCC", search_field="rownames", 
                           use_all_exp=FALSE, exp_name=NULL, case_insensitive=TRUE){
  if(!inherits(sce, "SummarizedExperiment")){
    stop("Only SummarizedExperiment object is supported (e.g., SummarizedExperiment, DESeqDataSet, or SingleCellExperiment)")
  }
  # Select which experiment to use
  if(use_all_exp){
    s <- unsplitAltExps(sce, prefix.rows=FALSE)
  }else{
    if(!is.null(exp_name)){
      if(exp_name %in% altExpNames(sce)){
        s <- altExp(sce, exp_name)
      }else{
        stop("Couldn't find the experiment name: '", exp_name, "'")
      }
    }else{
      s <- sce # Using main experiment
    }
  }
  
  # Validate search_field
  use_default_search_field = FALSE
  if(length(search_field) > 1){
    warning("Expecting the search_field to have length of 1, got ", length(search_field), " instead. The search_field will be set back to default setting (rownames(sce)).")
    use_default_search_field <- TRUE
  }
  
  if(search_field != "rownames"){
    # i.e., user specified colnames of rowData, checking if it's a valid choice
    if(ncol(rowData(sce)) == 0){
      warning("Empty rowData. The search_field will be set back to default setting (rownames(sce)).")
      use_default_search_field <- TRUE
      
    }else if(!(search_field %in% colnames(rowData(sce)))){
      warning("The given search_field doesn't match any column names of rowData(sce). The search_field will be set back to default setting (rownames(sce)).")
      use_default_search_field <- TRUE
    }
  }
  
  if(use_default_search_field){
    search_field="rownames"
  }
  
  
  # Search for the gene
  if(search_field=="rownames"){
    ref_vec <- rownames(s)
  }else{
    ref_vec <- rowData(s)[ , search_field]
  }
  
  if(case_insensitive){
    ref_vec <- tolower(ref_vec)
    query_regex <- tolower(query_regex)
  }
  
  hit_flag <- str_detect(ref_vec, query_regex)
  
  if(!any(hit_flag)){stop("Couldn't detect any query genes in the input")}
  
  s_fil <- s[hit_flag, ]
  return(s_fil)
}

sce_extract <- function(sce, query_regex="^ERCC", search_field="rownames", assay_names="counts", colData_colnames=character(0), 
                        use_all_exp=FALSE, exp_name=NULL, normalized=FALSE, case_insensitive=TRUE){
  # A wrapper function to run both sce_query_gene (filter genes) and sce_2_expr_tidy (filter assayNames, colData) in tandem
  # Filter SCE object for certain genes
  s_fil <- sce_query_gene(sce=sce, 
                          query_regex=query_regex, 
                          search_field=search_field, 
                          use_all_exp=use_all_exp, 
                          exp_name=exp_name, 
                          case_insensitive=case_insensitive)
  
  # Extract expression value
  expr_df <- sce_2_expr_tidy(s_fil, 
                             assay_names=assay_names, 
                             colData_colnames=colData_colnames, 
                             normalized=normalized)
  return(expr_df)
}


