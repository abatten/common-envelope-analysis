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
