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

### 5.plot_from_file.py
This script is for ploting the data created by the other analysis scripts.
It has options for plotting the energies, seperation, massloss and thermal.
To run the script you have to add speficy which plot you want by adding it as
an argument. 

#### Plot Types:
--energy: Create an energy plot
--seperation: Create a seperation plot
--angularmomentum: Create an angular momentum plot
--thermal: Create a thermal energy plot
--massloss: Create a mass loss and mass loss rate plot

#### Optionals
--smoothed: Find NaN's in energy file and replace them.
--marked: Add vertical black lines to the plot. List numbers after --marked. 

An example of this script being used to create a smooth, marked energy plot.
```
    python 5.plot_from_file.py /path/to/file/ --energy --smoothed --marked 1 2 3
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
