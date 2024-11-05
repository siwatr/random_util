# Transform fimo output to bed format
fimo2bed6(){
    # __NB__: Don't use this version of scripts
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
    fimo_file=$1
    min_pval=0.05
    min_qval=0.05
    
    # Reformat to BED6, then filter by p-value and q-value
    awk -v min_pval=$min_pval -v min_qval=$min_qval \
    'BEGIN{OFS="\t"}{ if ($8 <= min_pval && $9 <= min_qval ) print $3,$4,$5,$1"_"$2,$7,$6}' $fimo_file | grep -E "^chr" | sort -k1,1 -k2,2n > ${fimo_file%.tsv}.bed
    # awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$1"_"$2,$7,$6}' $fimo_file | grep -E "^chr" | sort -k1,1 -k2,2n > ${fimo_file%.tsv}.bed
    # # filter by p-value and q-value
    # awk -v min_pval=$min_pval -v min_qval=$min_qval '$8 <= min_pval && $9 <= min_qval' ${fimo_file%.tsv}.bed > ${fimo_file%.tsv}_filtered.bed
}

# awk -v min_pval=$min_pval -v min_qval=$min_qval 'BEGIN{OFS="\t"}{ if ( $8 >= min_pval && $9 >= min_qval ) print $3,$4,$5,$1"_"$2,$7,$6}' $fimo_file | grep -E "^chr" | sort -k1,1 -k2,2n | wc -l

# # Test
wd="/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/1_primary_analysis/motifsearch/JASPAR2022vertebrate_mm10"
# fimo2bed6 $wd/fimo/MA0062.3_GABPA/fimo.tsv
# cp $wd/fimo/MA0062.3_GABPA/fimo.bed $wd/../MA0062.3_GABPA_mm10.bed

# # Exclude empty lines and lines starting with "#"
# awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$1"_"$2,$7,$6}' $fimo_file | grep -E "^chr"





# Run from python script instead ==================================================================
python $wd/scripts/fimo_2_bed.py -h
# usage: fimo_2_bed.py [-h] [-i INPUT] [-o OUTFILE] [--no_filter] [-p MIN_PVAL] [-q MIN_QVAL] [-v]
# Convert FIMO output to BED format
# optional arguments:
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         FIMO output file
#   -o OUTFILE, --outfile OUTFILE
#                         Output directory (default: same directory as fimo_file)
#   --no_filter           Do not filter by p-value and q-value
#   -p MIN_PVAL, --min_pval MIN_PVAL
#                         Minimum p-value threshold (default: 0.05)
#   -q MIN_QVAL, --min_qval MIN_QVAL
#                         Minimum q-value threshold (default: 0.05)
#   -v, --verbose         Verbosity

python $wd/scripts/fimo_2_bed.py -i $wd/fimo/MA0062.3_GABPA/fimo.tsv -o $wd/bed/MA0062.3_GABPA_mm10.bed --no_filter

# list all directory that contain fimo.tsv in $wd/fimo
for f in $wd/fimo/*/fimo.tsv; do
    tf=$( basename $( dirname $f ) )
    # replace "-" with "_", then remove any duplicated "_"
    tf=$( echo $tf | tr "-" "_" | tr -s "_" )
    echo $tf
    python $wd/scripts/fimo_2_bed.py -i $f -o $wd/bed/${tf}_mm10.bed --no_filter
done

