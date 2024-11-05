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

gtf_path <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/0_refData/genome/mus_musculus/mm10/annotated_genes/ensembl/GRCm38.102/Mus_musculus.GRCm38.102_withChr.gtf"
gtf_path <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/0_refData/genome/mus_musculus/mm10/annotated_genes/custom_anno/C88delta_Zp3Cre/GRCm38.102_ERCC92_C88Delta.gtf"
if(!("gtf" %in% ls())){
  gtf <- BiocIO::import(gtf_path)
}
# Remove annoying GTF fields
mcols(gtf) <- as_tibble(mcols(gtf)) %>% 
  dplyr::select(type, gene_id, gene_name, transcript_id, transcript_name, exon_id)
gtf_gene <- gtf[gtf$type=="gene"]
strand_map <- dplyr::pull(as_tibble(gtf_gene), strand, name=gene_id)

wd <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/rna_seq/AI/P829_20240712_Nr5a2_cKO_polyA"
jcount_path <- paste0(wd, "/data/read_quant/star_cKO_featureCount_gene_body.txt.jcounts")
jcount_path <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/rna_seq/SR/SRR005_JG_C88Zp3Cre_GV/data/read_quant/star_multi_10_C88Delta_featureCount_gene_body.txt.jcounts"
sample_regex <- "(_Aligned.+bam)|(^\\.\\/)|(^\\.+)" # For renaming sample

jcnt <- read_featureCounts_jcounts(jcount_path, sample_regex=sample_regex, gene2strand_map=gtf_gene, rename_standard_col=TRUE)
jcnt_gr <- jcnt %>% 
  mutate(idx = seq_along(gene_id)) %>% 
  dplyr::filter(s1_seqnames == s2_seqnames) %>% 
  mutate(seqnames=s1_seqnames) %>% 
  # dplyr::select(-any_of(c("secondary_gene_id", "s1_seqnames", "s1_strand", "s2_seqnames", "s2_strand"))) %>% 
  dplyr::select(seqnames, start, end, strand, idx, gene_id) %>% 
  makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>% 
  suppressWarnings()


nr5a2_id <- gtf_gene[tolower(gtf_gene$gene_name)=="nr5a2"]$gene_id
jcnt_gr_nr5a2 <- jcnt_gr %>% 
  subset(., !is.na(.$gene_id)) %>% 
  subset(str_detect(.$gene_id, nr5a2_id))

jcnt_gr_nr5a2 <- junction_exon_map(junc_df=jcnt_gr_nr5a2, exon_gtf=gtf, return_GRanges=TRUE)
jcnt_gr_nr5a2
jcnt_gr_nr5a2_fil <- jcnt_gr_nr5a2[!is.na(jcnt_gr_nr5a2$start_exon) & !is.na(jcnt_gr_nr5a2$end_exon)]


# find where junction of interest are
# TODO: consider changing junction of interest to data frame instead
target_junc_delta <- list(left=c("ENSMUSE00000158720", "ENSMUSE00001343142"), right=c("ENSMUSE00000158724"))
target_junc_wt1 <- list(left=c("ENSMUSE00000158720", "ENSMUSE00001343142"), right=c("ENSMUSE00000158731"))
target_junc_wt2 <- list(left=c("ENSMUSE00000158731"), right=c("ENSMUSE00000158724"))

delta_idx <- jcnt_gr_nr5a2_fil[find_junction(jcnt_gr_nr5a2_fil, target_junc_delta)]$idx
wt1_idx <- jcnt_gr_nr5a2_fil[find_junction(jcnt_gr_nr5a2_fil, target_junc_wt1)]$idx
wt2_idx <- jcnt_gr_nr5a2_fil[find_junction(jcnt_gr_nr5a2_fil, target_junc_wt2)]$idx

int_idx <- c(delta_idx, wt1_idx, wt2_idx)
jcnt_nr5a2 <- jcnt[int_idx, ]
jcnt_nr5a2$junc_type = c("delta", "wt1", "wt2")
x <- jcnt_nr5a2[ ,  str_detect(colnames(jcnt_nr5a2), "(^S\\d+_)|(^junc_type$)")] %>% 
  column_to_rownames(var="junc_type") %>% 
  as.matrix() 

jcnt_nr5a2_tidy <-  jcnt_nr5a2 %>% 
  dplyr::select(gene_id, junc_type, dplyr::matches("S\\d+_Cre")) %>% 
  gather(key="sample", value="counts", -gene_id, -junc_type) %>% 
  mutate(condition=str_extract(sample, "(pos)|(neg)"))

jcnt_nr5a2_tidy %>% 
  ggplot(aes(x=condition, y=counts+1, color=condition)) +
  geom_boxplot(alpha=0.4) +
  # geom_violin(alpha=0.4) +
  geom_quasirandom(cex=2, alpha=0.5) +
  scale_y_log10() +
  facet_grid(.~junc_type) +
  theme_bw() +
  theme(legend.position="none")


x %>% 
  log10() %>% 
  replace(., is.infinite(.), 0) %>% 
  Heatmap(col=brewer.pal(n=8, "OrRd"), 
          cluster_columns=FALSE,
          column_split=str_extract(colnames(.), pattern="(neg)|(pos)"),
          top_annotation=HeatmapAnnotation(GT = str_extract(colnames(.), pattern="(neg)|(pos)")))




# tmp ----
jcnt_gr <- jcnt %>% 
  dplyr::filter(Site1_chr == Site2_chr) %>% 
  dplyr::select(seqnames=Site1_chr, start=Site1_location, end=Site2_location, name=PrimaryGene, strand) %>% 
  makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>% 
  suppressWarnings() 

# For IGV
export.bed(gr, con=paste0(wd, "/data/read_quant/star_cKO_featureCount_gene_body_spliceJunction.bed"))


# junction at nr5a2 genes

nr5a2_jcnt_gr <- jcnt_gr %>% 
  subset(!is.na(.$name)) %>% 
  subset(., str_detect(.$name, nr5a2_id))
nr5a2_jcnt_gr$idx <- seq_along(nr5a2_jcnt_gr)

tmp_gr <- nr5a2_jcnt_gr
tmp_gr$name <- seq_along(tmp_gr)
export.bed(tmp_gr, con=paste0(wd, "/data/read_quant/star_cKO_featureCount_gene_body_spliceJunction_NR5A2.bed"))

gtf_nr5a2_exon <- gtf %>% 
  subset(., .$type == "exon") %>% 
  subset(., seqnames(.)=="chr1") %>% 
  subset(., .$gene_name=="Nr5a2")
  # subset(., str_detect(tolower(.$gene_name), "nr5a2"))

# Remove transcript related field that make the data redundant
gtf_nr5a2_exon <- gtf_nr5a2_exon %>% 
  as_tibble() %>% 
  dplyr::select(seqnames, start, end, strand, any_of(c("type", "gene_id", "gene_name", "exon_id"))) %>% 
  unique() %>% 
  makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>% 
  BiocGenerics::sort()
gtf_nr5a2_exon$name <- gtf_nr5a2_exon$exon_id
# for IGV
export.bed(gtf_nr5a2_exon, con=paste0(wd, "/data/read_quant/tmp_nr5a2_exons.bed"))

# Find start exon overlap
start_ovl <- findOverlaps(flank(nr5a2_jcnt_gr, width=1, start=TRUE, both=TRUE), gtf_nr5a2_exon)
start_ovl_df <- as_tibble(start_ovl) %>% 
  dplyr::rename(q_hits = queryHits, s_hits = subjectHits) %>% 
  mutate(exon_id = gtf_nr5a2_exon$exon_id[s_hits]) %>% 
  group_by(q_hits) %>% 
  reframe(exon_id = paste0(unique(exon_id), collapse=";"))
nr5a2_jcnt_gr$start_exon <- NA
nr5a2_jcnt_gr$start_exon[start_ovl_df$q_hits] <- start_ovl_df$exon_id

# end exon overlap
end_ovl <- findOverlaps(flank(nr5a2_jcnt_gr, width=1, start=FALSE, both=TRUE), gtf_nr5a2_exon)
end_ovl_df <- as_tibble(end_ovl) %>% 
  dplyr::rename(q_hits = queryHits, s_hits = subjectHits) %>% 
  mutate(exon_id = gtf_nr5a2_exon$exon_id[s_hits]) %>% 
  group_by(q_hits) %>% 
  reframe(exon_id = paste0(unique(exon_id), collapse=";"))
nr5a2_jcnt_gr$end_exon <- NA
nr5a2_jcnt_gr$end_exon[end_ovl_df$q_hits] <- end_ovl_df$exon_id

nr5a2_jcnt_gr[!is.na(nr5a2_jcnt_gr$start_exon) | !is.na(nr5a2_jcnt_gr$end_exon)] %>% 
  as.data.frame() %>% 
  View()

# Check splicing at target junction
x <- gtf_nr5a2_exon %>% 
  as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, strand, type, gene_id, gene_name, transcript_id, transcript_name, exon_id) %>% 
  {split(., .$gene_name)}

unique(x$Nr5a2$exon_id)[!(unique(x$Nr5a2$exon_id) %in% unique(x$Nr5a2_Delta$exon_id))]

# Which exon does it connect to?

nr5a2_jcnt <- dplyr::filter(jcnt, str_detect(PrimaryGene, nr5a2_id))


# FUNCTIONS =======================================================================================
read_featureCounts_jcounts <- function(jcount_path, sample_regex=NULL, gene2strand_map=NULL, rename_standard_col=FALSE){
  jcnt <- as_tibble(read.table(jcount_path, header=TRUE, stringsAsFactors=FALSE))
  
  # Rename sample columns
  if(!is.null(sample_regex)){
    jcnt <- dplyr::rename_all(jcnt, .funs=str_remove_all, pattern=sample_regex)
  }
  
  if(!is.null(gene2strand_map)){
    if(inherits(gene2strand_map, "GRanges")){
      # Assume that GTF is provided (GRanges object)
      gtf <- gene2strand_map
      gene2strand_map <- as.character(strand(gtf))
      names(gene2strand_map) <- gtf$gene_name
    }
    
    jcnt$strand <- strand_map[jcnt$PrimaryGene]
  }
  
  # Re format standard column
  if(rename_standard_col){
    stdName_map <- c(PrimaryGene="gene_id",
                     SecondaryGenes="secondary_gene_id",
                     Site1_chr="s1_seqnames",
                     Site1_location="start",
                     Site1_strand="s1_strand",
                     Site2_chr="s2_seqnames",
                     Site2_location="end",
                     Site2_strand="s2_strand")
    colnames(jcnt)[colnames(jcnt) %in% names(stdName_map)] <- stdName_map[colnames(jcnt)[colnames(jcnt) %in% names(stdName_map)]]
  }
  
  return(jcnt)
}


junction_exon_map <- function(junc_df, exon_gtf, offset_bp=2, return_GRanges=FALSE){
  # Check inputs
  # junc_df: Junction counts data frame
  if(is.data.frame(exon_gtf)){
    junc_gr <- makeGRangesFromDataFrame(junc_df, keep.extra.columns=TRUE)
    
  }else if(inherits(exon_gtf, "GRanges")){
    junc_gr <- junc_df
  }
  
  # exon_gtf: GRanges object with at least one metadata column: exon_id
  required_cols <- c("type", "exon_id")
  err_msg <- paste0("Not all of the required columns: (", paste0(required_cols, collapse=", "), ") is present in the input.")
  if(is.data.frame(exon_gtf)){
    if(!all(required_cols %in% colnames(exon_gtf))){
      stop(err_msg)
    }
    exon_gtf <- makeGRangesFromDataFrame(exon_gtf, keep.extra.columns=TRUE)
    
  }else if(inherits(exon_gtf, "GRanges")){
    if(!all(required_cols %in% colnames(as.data.frame(mcols(exon_gtf))))){
      stop(err_msg)
    }
  }
  
  exon_gtf <- subset(exon_gtf, exon_gtf$type=="exon")
  mcols(exon_gtf) <- as_tibble(mcols(exon_gtf)) %>% 
    dplyr::select(any_of(required_cols))
  exon_gtf <- unique(exon_gtf)
  
  # overlap ---------------------------------------------------------------------------------------
  tmpFn_find_exon <- function(q, s, q_side="start", offset_bp=2){
    # Checking overlap between the ends side of the two features
    # i.e., comparing start to the end of the features
    q_gr <- flank(q, width=offset_bp, start=(q_side=="start"), both=TRUE)
    s_gr <- flank(s, width=offset_bp, start=(q_side!="start"), both=TRUE)
    ovl <- findOverlaps(q_gr, s_gr)
    
    # Capture which exon it being overlap with
    if(length(ovl)>0){
      ovl_df <- as_tibble(ovl) %>%
        dplyr::rename(q_hits = queryHits, s_hits = subjectHits) %>%
        mutate(exon_id = s_gr$exon_id[s_hits]) %>%
        group_by(q_hits) %>%
        reframe(exon_id = paste0(unique(exon_id), collapse=";"))
    }else{
      # Make empty output
      ovl_df <- tibble(qhits=0, exon_id="")
      ovl_df <- ovl_df[-1, ]
    }
    
    return(ovl_df)
  }
  
  start_exons <- tmpFn_find_exon(junc_gr, exon_gtf, q_side="start", offset_bp=offset_bp)
  end_exons <- tmpFn_find_exon(junc_gr, exon_gtf, q_side="end", offset_bp=offset_bp)
  junc_gr$start_exon <- NA
  junc_gr$end_exon <- NA
  if(nrow(start_exons)>0) {junc_gr$start_exon[start_exons$q_hits] <- start_exons$exon_id}
  if(nrow(end_exons)>0) {junc_gr$end_exon[end_exons$q_hits] <- end_exons$exon_id}
  
  if(return_GRanges){
    return(junc_gr)
  }else{
    return(as_tibble(junc_gr))
  }
}

find_junction <- function(junc_gr, target_junc){
  ## Example
  # target_junc <- list(left=c("ENSMUSE00000158720", "ENSMUSE00001343142"), 
  #                     right=c("ENSMUSE00000158724"))
  junc_list <- list(start=str_split(junc_gr$start_exon, ";"),
                    end=str_split(junc_gr$end_exon, ";"))
  
  # TODO: Consider using str_detect instead
  check1 <- sapply(junc_list$end, FUN=function(x){any(x %in% target_junc$left)})
  check2 <- sapply(junc_list$start, FUN=function(x){any(x %in% target_junc$right)})
  
  return(check1 & check2)
}

