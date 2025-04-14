#!/bin/bash

# List of IDs
ids=("A28600" "A28613" "A34434" "A26650" "A26675" "A28613" "A28641" "A28642" "A33846" "A34434" 
     "A26655" "A26659" "A26660" "A26671" "A28598" "A28607" "A28631" "A28636" "A29601" "A29611" 
     "A33838" "A33839" "A33844" "A33845" "A33847" "A33869" "A33880" "A33901" "A33904" "A33909" 
     "A33910" "A33911" "A33915" "A33925" "A34419" "A34423" "A34424" "A34429" "A34432")

ids_2=("A28600" "A28613" "A34434")

ids_204=("A34434"
"A26677"
"A29602"
"A33843"
"A33905"
"A34431"
"A34425"
"A34435"
"A28600"
"A33846")


# Create folders for each ID and subfolders AA197 to AA200
for id in "${ids[@]}"
do
    # Create the main folder for the ID
    mkdir "$id"
    
    # Create subfolders AA197 to AA200 within each ID folder
    for subfolder in 197-200
    do
        mkdir "$id/AA$subfolder"
    done
done


for id in "${ids_2[@]}"
do
    # Create the main folder for the ID
    mkdir "$id"

    # Create subfolders AA197 to AA200 within each ID folder
    for subfolder in 170-172
    do
        mkdir "$id/AA$subfolder"
    done
done

for id in "${ids_204[@]}"
do
    # Create the main folder for the ID
    mkdir "$id"

    # Create subfolders AA197 to AA200 within each ID folder
    for subfolder in 204
    do
        mkdir "$id/AA$subfolder"
    done
done


echo "Folders and subfolders created successfully!"
