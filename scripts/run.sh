#!/bin/bash

# The first argument is the script you want to run
SCRIPT=$1
INLIST=inlist_ce_analysis_template.ini

# If the script is the energy script
if [[ $SCRIPT == *"energy"* ]]; then
    echo " "
    echo "<------------>"
    echo "Running Energy Script"
    echo "<------------>"
    echo " "
    OUTPUT_NAMES=$"energy_components_"
    echo "Output File Names:" $OUTPUT_NAMES
fi

# Range depends on number of directories in path.
FINAL_PATH=5
for i in $(seq 0 ${FINAL_PATH})
do
    echo $i
    # Make Temp file to use as an inlist
    ANALYSIS_FILE=$(mktemp ce_analysis.XXXXX)
    # Read the inlist file and change the variables
    cat $INLIST |
    sed "s/INITIAL_PATH/$i/" |
    sed "s/FINAL_PATH_PLUS_ONE/$(($i+1))/" |
    sed "s/OUTPUT_FILE_NAME/$i/" > "$ANALYSIS_FILE"

    # Run the script for each directory at the same time.
    # /dev/null 2>&1 causes no output to be displayed.
    python $SCRIPT "$ANALYSIS_FILE" > /dev/null 2>&1 &
    rm $ANALYSIS_FILE
done
wait


