# Description: 
#    Read certain column of table (CSV or TSV) into an bash array variable.

# Dependencies:
#    load_flag_args.sh
#    - load_flag_args_boolean
#    - load_flag_args

# Usage:
#    read_table_to_array -A <array_name> -i <table_file> -sep <separator> -h -key <key_col> -val <val_col>
#        <array_name>: The name of the array variable to be created
#        <table_file>: The path to the table file
#        <separator>: The separator used in the table file. Default is tab*
#            * The conventional way of using $'\t' for tab separater won't work when it got parsed to load_flag_args. 
#              Thus we'll accept -sep '\t' instead, then convert it to a proper tab
#        <key_col>: The column index of the key. Default is 0, i.e., the first column
#        <val_col>: The column index of the value. If not set, it will be an indexed array. Otherwise, it will be an associative array.

# Read certain column of table (CSV or TSV) into an bash array variable.
# If val_col is not set, it will be an indexed array, otherwise it will be an associative array.
read_table_to_array(){
    # Load input arguments ------------------------------------------------------------------------
    # sep=$1; shift
    local array_name=$(load_flag_args "-A" "" $@)
    local table_file=$(load_flag_args "-i" "" $@)
    # The conventional way of using $'\t' for tab separater won't work when it got parsed to load_flag_args
    local sep=$(load_flag_args "-sep" "\t" $@)
    local header=$(load_flag_args_boolean "-h" "-header" $@)
    local key_col=$(load_flag_args "-key" "0" $@)
    local val_col=$(load_flag_args "-val" "" $@)
    
    # Check required input ------------------------------------------------------------------------
    # Check if array name is set
    if [[ -z "$array_name" ]]; then
        echo "Array name is not set"
        return 1
    fi
    
    # Check if the table file exists
    if [[ ! -s "$table_file" ]]; then
        echo "Table file does not exist"
        return 1
    fi
    
    # For debugging
    # echo "Key_col: $key_col | Val_col: $val_col"
    
    # Check if the separator is a tab
    if [[ $sep == "\t" || -z $sep ]]; then
        # echo "Tab separator"
    #    sep=$(echo $'\t')
       sep=$'\t'
    #    sep='    '
    fi
    
    # Check if val_col is not set
    if [ -z "$val_col" ]; then
        # echo "Value column is not set" # for debugging
        array_flag="a" # indexed array
    else
        array_flag="A" # associative array
    fi
    
    # a local function to print the table
    print_table(){
        local table_file=$1
        local header=${2:-"false"}
        if [ "$header" == "true" ]; then
            cat $table_file
        else
            tail -n +2 $table_file
        fi
    }
    
    # Build a recipe from ground up ---------------------------------------------------------------
    array_body=""
    while IFS=$sep read -r -a fields; do
        # Store value into an array. 
        key=${fields[$key_col]}
        # Behave differently depending on whether the value column is set or not
        if [ -z "$val_col" ]; then
            # Set into an indexed array
            array_body="${array_body} \"$key\""
        else
            value=${fields[$val_col]}
            # Store current key-value pair in the placeholder variable
            array_body="${array_body} [\"$key\"]=\"$value\""
        fi
    done < <(print_table "$table_file" $header) # Skip header row
    # Construct the entire recipe
    echo "declare -${array_flag} $array_name=( $array_body )"
}


# Example code:
if [ "${0}" = "${BASH_SOURCE}" ]; then
    # i.e., Script is being executed directly.
    # Example 1: no val is provided -- create an indexed array
    unset my_array_1 # make sure no array call my_array_1 exists
    read_table_to_array -A my_array_1 -i "./shell/read_table_to_array_example.tsv" -sep '\t' -key 0
    echo ${my_array_1[@]} # All elements in the array
    echo ${#my_array_1[@]} # Number of elements in the array
    echo ${!my_array_1[@]} # Name of the array

    # Example 2: val is provided -- create an associative array
    unset my_array_2
    read_table_to_array -A my_array_2 -i "./shell/read_table_to_array_example.tsv" -sep '\t' -key 0 -val 1
    echo ${my_array[@]}
    echo ${#my_array[@]}
    echo ${!my_array[@]}

    # Example 2.2: Reading CSV file
    unset my_array_2
    read_table_to_array -A my_array_2 -i "./shell/read_table_to_array_example.tsv" -sep ',' -key 0 -val 1
    echo ${my_array[@]}
    echo ${#my_array[@]}
    echo ${!my_array[@]}
    
else
    # echo "Script is being sourced."
fi

