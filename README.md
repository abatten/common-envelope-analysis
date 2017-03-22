# yt-scripts
My collection of yt-scripts that I have written for data analysis. 
These scripts have been specifically written to analyse the output 
of an ENZO common envelope simulation.



## Analysis Scripts:
### 1.damping_analysis.py
This script is designed to analyse the data of a pre-common
envelope damping simulation. With a single star that is inserted 
in enzo from a 1D stellar evolution code output, via the YT package.

The script can generate a number of plots as listed below. 
To run this script you need to specify an inlist file with all the 
required parameters and specify the plots by adding them as 
arguments when calling the python script. 

For example
```
    python 1.damping_analysis.py inlist_ce_analysis 0,9,18
```
### 2.ce_energy_smoothed_potential.py
```
    python 2.ce_energy_smoothed_potential.py inlist_ce_analysis
```

### 3.ce_angular_momentum.py
```
    python 3.ce_angular_momentum.py inlist_ce_analysis
```
### 4.ce_seperations.py
```
    python 4.ce_seperations.py inlist_ce_analysis
```
## CEFunctions
CEFunctions is a modules that is required to run these YT scripts. 
Each script will import CEFunctions as cef.
```
import cefunctions as cef
```
Within CEFunctions you will find the following functions:

**index2str:** Converts an index to a string with any required number of pre-fixes.
```
index2str(4) --> 0004
```
**root_sort:** Takes an input directory and sorts all the common envelope data files.

**distance:** Calculates the distance between two 3D vectors and also 
allows for unit conversion. 

**grav_pot:** Calculates the gravitational potential energy between two masses. 
grav_pot also has the functionality to use a smoothed gravitational potential energy 
as given by M.Ruffert 1993.
