#!/bin/bash

# This script is used to run multiple instances of the same
# Common Envelope python script to help speed up the analysis.

# The first argument is the script you want to run
SCRIPT=$1

# This is the inlist of parameters that will be modified by this script.
INLIST=inlist_ce_analysis_template.ini

# Range depends on number of directories in path.
# The FINAL_PATH is equal to the total number of directories.
FINAL_PATH=176

for i in $(seq 0 ${FINAL_PATH})
do
    echo $i
    # Create a temp inlist for the python script to read
    ANALYSIS_FILE=$(mktemp ce_analysis.XXXXX)
    # Read the inlist file and change the variables
    cat $INLIST |
    sed "s/INITIAL_PATH/$i/" |
    sed "s/FINAL_PATH_PLUS_ONE/$(($i+1))/" |
    sed "s/OUTPUT_FILE_NAME/$i/" > "$ANALYSIS_FILE"

    # Run the script for each directory at the same time.
    # /dev/null 2>&1 causes no output to be displayed.
    nice -n 19 python $SCRIPT "$ANALYSIS_FILE" > /dev/null 2>&1 &
done
wait

# Remove all the temp files after the script finishes.
rm ce_analysis.*

echo " "
echo "FINISHED CE_RUN: " $SCRIPT
echo " "
