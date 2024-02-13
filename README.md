# MC_Hopper
*by Jolla Kullgren*

A simple lattice kinetic Monte-Carlo program written in Fortran90. 

The source code is contained in a single file called `kMC.f90` and the main input files are `input.dat` and `geom.dat`.

## Geometry - Lattice specification 
The file `geom.dat` has the following format:

The first line specifies the number of points, and the following lines have the format:
```
x y z type_of_site type_of_species No_of_neigbors No_of_neigbors*id No_of_neigbors*jump_type
```
`No_of_neigbors` gives the number of neigboring sites including the site it self (some events are "selfevents"). `id` are the id number of the neigbor. Id numbers are defined from 1 to the number of sites as they are listed in the file. `jump_type` specifies the type of connection that the two sites has. This is useful for defining long and short migration steps etc.

## Barriers and general input
The file `input.dat` has the following format:

First line:
```
No_of_2body_interactions No_of_1body_interactions No_of_barriers No_of_simulation_steps Printing_freq printing_type Temperature
```
second line:
```
No_of_species_types No_of_site_types No_of_jump_types
```


the following `No_of_2body_interactions` lines:
```
jump_type site_A site_B  spec_at_A spec_at_B Energy
```

the following `No_of_1body_interactions` lines:
```
site  spec Energy
```

the following `No_of_barriers` lines:

```
jump_type site_A site_B init_spec_at_A init_spec_at_B final_spec_at_A final_spec_at_B Attempt_freq Barrier
```

The barrier for a general event involving two sites `site_A` and `site_B` connected with `jump_type`. The initial conditions to be fulfilled are specified by
`init_spec_at_A` and `init_spec_at_B` and the effect of the event by `final_spec_at_A` and `final_spec_at_B`.

## General output

The accumulated time and the lifetime of the current frame. These are printed in the standard output from the code. 


The evolution of the system is recorded in the file `KMC.xyz`. The first line of this file gives the the number of recorded frames after which blocks describing these frames follows. The blocks first give the number of lattice sites. Each lattice site, in each of the frames, are represented by an integer number specifying the species occupying the lattice site followed by the corresponding cartesian coordinate to that site.



