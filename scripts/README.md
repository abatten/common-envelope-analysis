## Analysis Scripts:
### ce_damping_analysis.py
This script is designed to analyse the data of a pre-common
envelope damping simulation. With a single star that is inserted 
in enzo from a 1D stellar evolution code output, via the YT package.

The script can generate a number of plots as listed below. 
To run this script you need to specify an inlist file with all the 
required parameters and specify the plots by adding them as 
arguments when calling the python script. 


#### List of Avaliable Plots:

0: Density vs Radius          
1: Radial Velocity vs Radius        
2: Kinetic Energy vs Radius   
3: Thermal Energy vs Radius   
4: Total Energy vs Radius     
5: Grav Potential vs Radius   
6: X-Velocity vs Radius        
7: Y-Velocity vs Radius       
8: Z-Velocity vs Radius       
9: Density Slice along Z-Axis                   
10: Density Projection along Z-Axis    
11: Pressure Slice along Z-Axis     
12: Thermal Energy SLice along Z-Axis       
13: Density Gradient Modulus Slice along Z-Axis   
14: Velocity Modulus Slice along Z-Axis      
15: Temperature Slice along Z-Axis     
16: Mach Number Slice along Z-Axis      
17: Entropy Slice along Z-Axis       
18: Gravitational Potential Slice along Z-Axis     

An example of running this script to create density and gravitational potential
slice plots.
```
    python ce_damping_analysis.py inlist_ce_analysis 0,9,18
```
### ce_energy_smoothed_potential.py
This script calculats all of the energy components of the common envelope
```
    python ce_energy_smoothed_potential.py inlist_ce_analysis.ini
```

### ce_angular_momentum.py
```
    python ce_angular_momentum.py inlist_ce_analysis.ini
```
### ce_seperations.py
This script calculates the seperation(s) between the primary and the companion(s).
```
    python ce_seperations.py inlist_ce_analysis.ini
```

### ce_plot_from_file.py
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
--smoothed: Find NaN's in energy file and replace them  
--marked: Add vertical black lines to the plot. List numbers after --marked   

An example of this script being used to create a smooth, marked energy plot.
```
    python ce_plot_from_file.py /path/to/file/ --energy --smoothed --marked 1 2 3
```
