# MC_Hopper
*Written by Jolla Kullgren*

Lattice Monte-Carlo and kinetic Monte-Carlo programs written in Fortran90. 
The main input files are `input.dat` and `geom.dat`.

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


