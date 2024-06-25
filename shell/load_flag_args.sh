# Description: 
#    A utility function to looking for a flag in the command line arguments, and return the value that comes after it.
# Usage:
#    load_flag_args_boolean <short_flag> <long_flag> $@
#        <short_flag>: The short flag to look for, e.g., "-h"
#        <long_flag>: The long flag to look for, e.g., "--help"

# Cautions:
#   - It cannot distinguish a collection of short flags (e.g., -vf (same as -v -f)) like a typical bash command.
#   - it just looking for a literal match of the flag, without really check if it's the flag or not, so one could also put anything in there.

# A utility function to load command line arguments of certain flag
load_flag_args_boolean(){
  # flag to looking for, e.g. "-h"
  local flag_short=$1; shift
  local flag_long=$1; shift
  
  # Find the flag and instore the value that comes after until the new flag is found
  flag_found=false
  while [[ $# -gt 0 ]]; do
    key="$1"
    # Check if the key is the flag we want to read
    if [[ $key == "$flag_short" || $key == "$flag_long" ]]; then
      flag_found=true
      break # no need to do anything else
    fi
    shift
  done
  # Return output
  echo $flag_found
}


# Usage:
#    load_flag_args <flag> <default_value> $@
load_flag_args(){
  local flag=$1; shift
  local default_value=$1; shift
  # echo $@ # for debugging
  
  # Find the flag and instore the value that comes after until the new flag is found
  declare -a bw_keys # a placeholder
  read_mode=false
  while [[ $# -gt 0 ]]; do
    key="$1"
    # Check if the key is the flag we want to read
    if [[ $key == "$flag" ]]; then
      read_mode=true
      shift
      continue # go the the next iteration
    elif [[ $key == "-"* ]]; then
      # Turn off the read mode
      read_mode=false
    fi

    if [[ $read_mode == true ]]; then
      bw_keys+=($key)
    fi
    shift
  done

  # Set default values:
  if [ ${#bw_keys[@]} -eq 0 ]; then
    bw_keys=($default_value)
  fi
  # Return output
  echo ${bw_keys[@]}
}


# Example code:
if [ false == true ]; then
    # Boolean parameter test
    arg_val="A B C"
    load_flag_args_boolean "-v" "--verbose" -a 1 -b -c x y z -d $arg_val # no -v flag
    load_flag_args_boolean "-v" "--verbose" -a 1 -b -c x y z -d $arg_val # -v in the comment shouldn't be counted
    load_flag_args_boolean "-v" "--verbose" -v -a 1 -b -c x y z -d $arg_val # -v at the beginning
    load_flag_args_boolean "-v" "--verbose" -a 1 -v -b -c x y z -d $arg_val # -v in the middle
    load_flag_args_boolean "-v" "--verbose" -a 1 -b -c x y z -d $arg_val -v # -v at the end
    load_flag_args_boolean "-v" "--verbose" -a 1 -b -c v x y z -d $arg_val # v appear in other parameter
    load_flag_args_boolean "-v" "--verbose" --verbose -a 1 -b -c x y z -d $arg_val # Full flag
    load_flag_args_boolean "-v" "--verbose" -verbose -a 1 -b -c x y z -d $arg_val  # Full flag with only one hyphen
    load_flag_args_boolean "-v" "--verbose" --verbose1 -a 1 -b -c x y z -d $arg_val # Mismatched flag
    # NB: it cannot distinguish a collection of short flags like a typical bash command.
    load_flag_args_boolean "-v" "--verbose" -va -a 1 -b -c x y z -d $arg_val # Mismatched flag

    # Regular parameter test
    load_flag_args "-c" "default" -a 1 -b -d $arg_val            # no -c flag
    load_flag_args "-c" "default" -a 1 -b -c x -d $arg_val       # -c flag with a single values
    load_flag_args "-c" "default" -a 1 -b -c x y z -d $arg_val   # -c flag with multiple values
    load_flag_args "-c" "default" -a 1 -b -c -d $arg_val         # -c flag without value
fi