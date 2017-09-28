# Common Envelope Analysis
These scripts are designed to analyse the data of a common
envelope simulation. With a single star that is inserted 
in ENZO from a 1D stellar evolution code output, via the YT package.

The simulations could have one or more companion particles and these
scripts account for the number of particles.

## Inlists:
These common envelope analysis scripts have been writen to require an
'inlist' to specify the individual parameters of the scripts.

`inlist_ce_analysis.ini`

The inlist is broken up into sections. The first section is a 
'Common Section' which every script will read such as the `root_dir` or
`plot_dir` which is are variables for where the data is stored and where 
to write output files respectfully.

After the common section are the remaining sections that are specific
to each script. For example; `ce_energy_smoothed_potential` will only read
the common section and the energy section. These parameters can be changed
as needed.

## Multirun
Due to large number of numerical calculations some of these scripts
perform (in particular `ce_energy_smoothed_potential.py` and 
`ce_mass_loss.py`) they can take a long time to run in sequence. The 
shell script `multirun.sh` is basically an easy parallelisation 
of my code. It works by launching multiple instances of the python script
on different sets of data simulataneiously. 

To launch a script using multirun first modify the multirun inlist:

```
inlist_ce_analysis_multirun.ini
```

Then launch `multirun` parsing the python script as the first argument

```
    ./multirun.sh ce_energy_damping_analysis.py
```
This will create an output file for every data dump. To merge the outputs
into a single file use `ce_file_merge.py`.

## Analysis Scripts:
|Script Name                    |Purpose                                            |
|---                            |---                                                |
|`ce_angular_momentum`          |Compute angular momentum components                |
|`ce_centre_of_mass`            |Find the centre of mass of the system              |
|`ce_core_seperations`          |Compute the seperation between particles           |
|`ce_damping_analysis`          |Produce various radial/slice plots from a dump     |
|`ce_energy_smoothed_potential` |Compute the energy components                      |
|`ce_file_merge`                |Merge files together after `multirun`              |
|`ce_gas_outflow`               |Find regions where gas is outflowing               |
|`ce_gravodrag`                 |Compute the gravodrag on the companions            |
|`ce_mass_loss`                 |Compute bound/unbound mass of the envelope         |
|`ce_orbits`                    |Find the orbit of the particles from velocity      |
|`ce_particle_pos_and_vel`      |Find the positions and velocities of the particles |
|`ce_plot_from_file`            |Read a text file of data and produce various plots |
|`ce_vel_vs_kep.py              |Compare particle velocities to keplerian velocity  |

### ce_angular_momentum.py
Calculates the angular momentum components of the companions.

Note: Still under work
```
    python ce_angular_momentum.py inlist_ce_analysis.ini
```
#### Inlist Parameters
output_file_name (str) : Name of the file

output_file_append (bool) : Append to file if existing

angular_momentum_wrt_com (bool) : Find Angular momentum with respect 
to centre of mass

### ce_centre_of_mass.py
Find the centre of mass in the box

### ce_core_seperations.py
Calculates the seperation(s) between the primary and the companion(s).

```
    python ce_seperations.py inlist_ce_analysis.ini
```

### ce_damping_analysis.py
Can generate a number of plots as listed below. 

|Radial Plots                   |Axial Plots                               |
|---                            |---                                       |
|0: Density vs Radius 		    |9: Density Slice Z-Axis   	               |
|1: Radial Velocity vs Radius   |10: Density Projection Z-Axis   	       |
|2: Kinetic Energy vs Radius	|11: Pressure Slice Z-Axis  	           |
|3: Thermal Energy vs Radius	|12: Thermal Energy Slice Z-Axis   	       |
|4: Total Energy vs Radius      |13: Density Gradient Modulus Slice Z-Axis |
|5: Grav Potential vs Radius    |14: Velocity Modulus Slice Z-Axis         |
|6: X-Velocity vs Radius        |15: Temperature Slice Z-Axis              |
|7: Y-Velocity vs Radius        |16: Mach Number Slice Z-Axis              |
|8: Z-Velocity vs Radius        |17: Entropy Slice Z-Axis                  |
|                               |18: Gravitational Potential Slice Z-Axis  |

To run this script you need to specify the plots by adding them as 
arguments when calling the python script. 

An example of running this script to create density and gravitational potential
slice plots.
```
    python ce_damping_analysis.py inlist_ce_analysis.ini 9,18
```

### ce_energy_smoothed_potential.py
Calculates all of the energy components of the common envelope. 
It accounds for ENZO using a smoothed potential at small seperations. It
also automatically excludeds cells with potential energy spikes due to AMR.

```
    python ce_energy_smoothed_potential.py inlist_ce_analysis.ini
```

### ce_mass_loss.py
Calculates the mass lost from the system through mass leaving the grid.

```
    python ce_mass_loss.py inlist_ce_analysis.ini
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
--posvel: Create a positions and velocities plot
--gravodrag: Create plot of Gravitational drag


#### Optionals
--smoothed: Find NaN's in energy file and replace them.  
--marked: Add vertical black lines to the plot. List numbers after --marked   

An example of this script being used to create a smooth, marked energy plot.
```
    python ce_plot_from_file.py /path/to/file/ --energy --smoothed --marked 1 2 3
```

### ce_vel_vs_kep.py
Calculates the velocities and seperations of the particles in the simulations 
and compares these velocities to the Kelperian velocity at that seperation.

Currently Under Work
