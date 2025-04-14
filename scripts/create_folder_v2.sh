#!/bin/bash

# List of libraries
libraries=(
    "Lib.K.2436"
    "Lib.K.2437"
    "Lib.K.2438"
    "Lib.K.2440"
    "Lib.K.2441"
    "Lib.K.2442"
    "Lib.K.2443"
    "Lib.K.2445"
    "Lib.M.2186"
    "Lib.M.2189"
    "Lib.M.2190"
    "Lib.M.2191"
    "Lib.M.2192"
    "Lib.M.2193"
    "Lib.M.2194"
    "Lib.M.2195"
    "Lib.M.2196"
    "Lib.M.2197"
)

# Loop over each library and create the directories
for lib in "${libraries[@]}"
do
    # Create the main library directory
    mkdir -p "$lib/TW1+TW2"
done

echo "Folders created successfully!"
