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
    
    # Print input arguments
    printf "DT_wrapper_2 Inputs:\n---------------------------\n"
    printf "\tprefix : $prefix\n"
    printf "\tforce_run : $force_run\n"
    printf "\tbw_keys: ${bw_keys[@]}\n"
    printf "\tbw_info: ${bw_info[@]}\n"
    printf "\tbed_file: ${bed_file[@]}\n"
    printf "\n"
    
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
    
    # Non-bw-associated arguments -----------------------------------------------------------------
    ## TODO: Consider turning this into a function
    
    # Detect additional arguments for computeMatrix
    # Set the default value for some of the arguments
    declare -A accept_comMat_args=( \
     [-p]=16 \
     ["--upstream"]=2000 ["--downstream"]=2000 ["--binSize"]=100 \
     ["--referencePoint"]="center" \
    )
    
    comMat_args=""
    for i in "${!accept_comMat_args[@]}"; do
        # echo "$i: ${accept_comMat_args[$i]}" # For debugging
        # Check if user provide such argument
        local a=$(load_flag_args $i "${accept_comMat_args[$i]}" $@)
        # They should all be non empty variable, but just in case:
        if [ ! -z "$a" ]; then
            comMat_args="$comMat_args $i $a"
        fi
    done
    
    # Check for boolean arguments
    # NB: don't activate for now, need to think about how to deal with the default value here
    # declare -A accept_comMat_args_boolean=( ["--missingDataAsZero"]="true" ["--skipZeros"]="true" )
    # for i in "${!accept_comMat_args_boolean[@]}"; do
    #     local a=$(load_flag_args_boolean $i ${accept_comMat_args_boolean[@]} $@)
    #     if [ "$a" == true ]; then
    #         comMat_args="$comMat_args $i"
    #     fi
    # done
    
    # Detect additional arguments for plotHeatmap
    # Set the default value for some of the arguments
    declare -A accept_PH_args=( \
     ["--kmeans"]=4 ["--clusterUsingSamples"]=1 \
     ["--refPointLabel"]=${accept_comMat_args["--referencePoint"]} \
    )

    PH_args=""
    for i in "${!accept_PH_args[@]}"; do
        # echo "$i: ${accept_PH_args[$i]}" # For debugging
        # Check if user provide such argument
        local a=$(load_flag_args $i "${accept_PH_args[$i]}" $@)
        # They should all be non empty variable, but just in case:
        if [ ! -z "$a" ]; then
            PH_args="$PH_args $i $a"
        fi
    done
    
    ## Load arguments associated with bw inputs ---------------------------------------------------
    ## TODO: Consider turning this into a function
    
    ## Search the column index for value from the header instead of explicitly providing it
    # bw_keys="Nr5a2 ATAC_GSE66390 H3K9me3_GSE97778_rep3" # Test param
    
    # Set column starting index based on table reader function that we'll use
    #   read_table_array() use native bash indexing (start with 0), 
    #   while read_table_array_2() use awk indexing (start with 1)
    local col_start_idx=1 # Where does the first column start?
    
    ## Read header of the bw_info file
    IFS=$'\t' read -r -a bw_header <<< $(head -n 1 $bw_info)
    # unset bw_args
    declare -a arg_name=( "label" "path" "zMax" "yMax" "yMin" "colorMap" )
    declare -A bw_args # The final arguments for each specified bw sample input
    for i in ${arg_name[@]}; do
        val_col_idx=""
        # find column index within the header
        for h in "${!bw_header[@]}"; do
            # echo "$h: ${bw_header[$h]}"
            if [[ ${bw_header[$h]} == $i ]]; then
                val_col_idx=$h
                break
            fi
        done
        
        # Check if the column is found
        if [ -z "$val_col_idx" ]; then
            echo "Column $i not found in the header of the $(basename $bw_info) file"
            bw_args[$i]="" # Just in case, set it to empty string for now
        else
            # Column found: load the values
            # unset "bw_${i}_array"
            eval $(read_table_to_array_2 -A "bw_${i}_array" -i $bw_info -sep '\t' -key 1 -val $((val_col_idx+$col_start_idx)) )
            bw_args[$i]=$(subset_bash_array "bw_${i}_array" $bw_keys)
            # echo "$i:    ${bw_args[$i]}" # debug
        fi
    done
    
    # Running DeepTools commands ------------------------------------------------------------------
    #### Compute matrix
    # Check if the matrix output file exists, if not, run the computeMatrix
    if [[ ! -s $DT_dir/matrix/${prefix}.gz || $force_run == true ]]; then
        # Assuming that user provide path and label in bw_args
        cmd="computeMatrix reference-point $comMat_args \
         --missingDataAsZero --skipZeros \
         -R $bed_file -S ${bw_args["path"]} --samplesLabel ${bw_args["label"]} \
         -o $DT_dir/matrix/${prefix}.gz"
        
        # Show the command and run
        echo $cmd
        printf "\n"
        eval $cmd
    else
        echo "Matrix file already exists, skipping computeMatrix step."
    fi
    
    #### Plotting heatmap
    # TODO: Consider turning this into a function
    # Validate bw-associated arguments required for plotHeatmap (PH)
    PH_bw_arg_name=("zMax" "yMax" "yMin" "colorMap")
    PH_bw_args=""
    for i in "${PH_bw_arg_name[@]}"; do
        cur_arg=${bw_args[$i]} # Extract input argument
        # If the string is not empty, add it to the command
        if [[ ! -z $cur_arg ]]; then
            PH_bw_args="${PH_bw_args} --${i} ${cur_arg}"
        fi
    done
    
    ph_cmd="plotHeatmap -m $DT_dir/matrix/${prefix}.gz -out $DT_dir/plots/${prefix}_heatmap.png $PH_args $PH_bw_args"
    
    # Show the command and run
    echo $ph_cmd
    printf "\n"
    eval $ph_cmd
}
