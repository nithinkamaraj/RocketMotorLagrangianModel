#!/bin/bash
:set ff=unix
# Determine the directory where the script is located
script_dir="${0%/*}"

# Go one level up
one_level_up="${script_dir%/*}"

# Print the directory
echo "The script is located in: $one_level_up"
