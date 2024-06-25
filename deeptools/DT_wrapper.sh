# Description:
#     A wrapperfor running deeptools computeMatrix and plotHeatmap

# Dependencies:
# load_flag_args.sh
#    - load_flag_args_boolean
#    - load_flag_args
# subset_bash_array.sh
# read_table_to_array.sh

DT_wrapper(){
    local prefix=$1
    shift
    local bed_file=$@
    # Todo: 
    # - add force run arg
    # - automatically read bw input from file
    echo $(basename $bed_file)
    # Will use most of its arguments from the global variables
    computeMatrix reference-point -p 16 \
     -R $bed_file \
     -S $bw_path_args \
     --samplesLabel $bw_labl_args \
     --upstream 2000 --downstream 2000 --binSize 50 --missingDataAsZero --skipZeros \
     -o $DT_dir/matrix/${prefix}.gz
    
     plotHeatmap -m $DT_dir/matrix/${prefix}.gz -out $DT_dir/plots/${prefix}_heatmap.png \
       --kmeans 4 --clusterUsingSamples 1 \
       --zMax $bw_zMax_args --yMin $bw_yMin_args --yMax $bw_yMax_args \
      --colorMap "Greys" --refPointLabel "B1 center"
}

# =================================================================================================
# A wrapperfor running deeptools computeMatrix and plotHeatmap
DT_wrapper_2(){
    # Get input arguments ---------------------------------------------------------------------------
    local prefix=$1; shift
    local force_run=$(load_flag_args_boolean "-f" "--force" $@)
    local bed_file=$(load_flag_args "-bed" "" $@)
    local bw_info=$(load_flag_args "-bw_info" "" $@)
    local bw_keys=$(load_flag_args "-bw" "" $@)
    local upstream=$(load_flag_args "--upstream" "2000" $@)
    local downstream=$(load_flag_args "--downstream" "2000" $@)
    local bin_size=$(load_flag_args "--binSize" "100" $@)
    
    printf "DT_wrapper_2 Inputs:\n---------------------------\n"
    printf "\tprefix : $prefix\n"
    printf "\tforce_run : $force_run\n"
    printf "\tbw_keys: ${bw_keys[@]}\n"
    printf "\tbw_info: ${bw_info[@]}\n"
    printf "\tbed_file: ${bed_file[@]}\n"
    printf "\n"
    # Todo: 
    # - add force run arg
    # Check required input ------------------------------------------------------------------------
    # Check if the prefix exist
    if [ -z "$prefix" ]; then
        echo "Prefix is not set"
        return 1
    fi
    
  # Check if the bw_info exists
    if [ ! -s "$bw_info" ]; then
        echo "bw_info file does not exist"
        return 1
    fi
    
    # Load arguments associated with bw inputs
    # bw_keys="Nr5a2 ATAC_GSE66390 H3K9me3_GSE97778_rep3"
    eval $(read_table_to_array -A bw_file_array -i $bw_info -sep '\t' -key 0 -val 1)
    eval $(read_table_to_array -A bw_zMax_array -i $bw_info -sep '\t' -key 0 -val 2)
    eval $(read_table_to_array -A bw_yMax_array -i $bw_info -sep '\t' -key 0 -val 3)
    eval $(read_table_to_array -A bw_yMin_array -i $bw_info -sep '\t' -key 0 -val 4)
    eval $(read_table_to_array -A bw_colMap_array -i $bw_info -sep '\t' -key 0 -val 5)
    
    local bw_labl_args=$bw_keys
    local bw_path_args=$(subset_bash_array bw_file_array $bw_keys)
    local bw_zMax_args=$(subset_bash_array bw_zMax_array $bw_keys)
    local bw_yMax_args=$(subset_bash_array bw_yMax_array $bw_keys)
    local bw_yMin_args=$(subset_bash_array bw_yMin_array $bw_keys)
    local bw_colMap_args=$(subset_bash_array bw_colMap_array $bw_keys)
    
    # For debugging
    # echo $bw_path_args
    # echo $(basename $bed_file)
    
    # # Will use most of its arguments from the global variables
    # Check if the matrix output file exists, if not, run the computeMatrix
    if [[ ! -s $DT_dir/matrix/${prefix}.gz || $force_run == true ]]; then
        cmd="computeMatrix reference-point -p 16 \
         --referencePoint center \
         -R $bed_file \
         -S $bw_path_args \
         --samplesLabel $bw_labl_args \
         --upstream $upstream --downstream $downstream --binSize $bin_size --missingDataAsZero --skipZeros \
         -o $DT_dir/matrix/${prefix}.gz"
        
        # Print cmd and run
        echo $cmd
        eval $cmd
    else
        echo "Matrix file already exists, skipping computeMatrix step."
    fi
    
    cmd="plotHeatmap -m $DT_dir/matrix/${prefix}.gz -out $DT_dir/plots/${prefix}_heatmap.png \
     --kmeans 4 --clusterUsingSamples 1 \
     --zMax $bw_zMax_args --yMin $bw_yMin_args --yMax $bw_yMax_args \
     --colorMap $bw_colMap_args --refPointLabel B1 center"
    echo $cmd
    eval $cmd
}
