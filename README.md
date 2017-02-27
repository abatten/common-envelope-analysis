# yt-scripts
My collection of yt-scripts that I have written for data analysis.

## Analysis Scripts:
1.damping_analysis

2.ce_energy_smoothed_potential

3.ce_angular_momentum

## CEFunctions
CEFunctions is a modules that is required to run these YT scripts. 
You can import CEFunctions as follows. 
--> import cefunctions as cef

Within CEFunctions you will find the following functions:

index2str: Converts an index to a string with any required number of pre-fixes.
		e.g index2str(4) --> 0004


root_sort: Takes an input directory and sorts all the common envelope data files.
