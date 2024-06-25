# Description: This script is used to subset an associative array in bash
# Usage: subset_bash_array <array_name> <keys>
#    <array_name>: An associative array available in the current scope/environment
#                  Only the name of the variable is required -- i.e., no need to use $ or {}
#   <keys>: A list of keys to be selected from the associative array
subset_bash_array(){
    ref_array=$1
    # Select elements from the associative array
    # Use nameref to refer to the global associative array, allowing us to access the element of ref_array
    #     shift: Remove the first argument of this function, so that the rest are the keys
    declare -n ref_array=$1; shift  
    selected_keys=("$@")
    
    local val=""
    for key in "${selected_keys[@]}"; do
        val="$val ${ref_array[$key]}"
        # echo "${ref_array[$key]} "
    done
    
  echo ${val[@]}
}

# Simialar to the previous one, but return the script to generate a new array instead of printing the values
subset_bash_array_2(){
    # ref_array=$1
    declare -n ref_array=$1; shift  
    new_array=$1; shift
    selected_keys=("$@")
    
    local array_body=""
    for key in "${selected_keys[@]}"; do
        array_body="$array_body [\"$key\"]=\"${ref_array[$key]}\""
        # echo "${ref_array[$key]} "
    done
    
    # Making an associative array as a default behavior (at least for now)
    array_cmd="declare -A $new_array=( $array_body )"
    echo $array_cmd
}


# Example code: 
if false=true; then
    # Define an associative array
    declare -A my_array
    my_array[key1]="value1"
    my_array[key2]="value2"
    my_array[key3]="value3"
    
    # Call the function
    subset_bash_array my_array "key1" "key3"
    # Output: value1 value3
    subset_bash_array_2 my_array new_array "key1" "key3"
    # declare -A new_array=( ["key1"]="value1" ["key3"]="value3" )
    
    # we'll need to eval it to assign to the new array
    unset new_array
    eval $(subset_bash_array_2 my_array new_array "key1" "key3")
    echo ${!new_array[@]}
    echo ${new_array[@]}
fi

# DEPRECIATED
# NB: "${0}" = "${BASH_SOURCE}" is used to check whether a script is being executed directly or being sourced from another script or the command line
# ${0} is the name of the script itself
# ${BASH_SOURCE} is the name of the script that is being sourced

