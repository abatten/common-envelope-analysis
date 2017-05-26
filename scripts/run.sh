#!/bin/bash


# The first argument is the script you want to run
SCRIPT=$1

# Range depends on number of directories in path.
for i in {0..142}
do
    echo $i
    # Read the inlist file and change the variables
    cat inlist_ce_analysis_template2.txt |
    sed "s/INITIAL_PATH/$i/" |
    sed "s/FINAL_PATH_PLUS_ONE/$(($i+1))/" |
    sed "s/OUTPUT_FILE_NAME/$i/" > inlist_ce_analysis_template2.txt

    # Run the script for each directory at the same time.
    # /dev/null 2>&1 causes no output to be displayed.
    python $SCRIPT inlist_ce_analysis_template2.txt & #> /dev/null 2>&1 &
done
wait
