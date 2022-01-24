These scenarios were generated using poets-dpd

Scenarios
---------

### Water

This is just plain water, so no bonds, all beads interact in the same way.
There are two bead-types, just for visualisation purposes, but it is all
water. Volume is scaled as a cube.

### Water_density

Plain cubes of water, but the density is varied from 0.5 through to 8.0
in increments of 0.5. These are available in cubes of side 25, 35, 45, and 55.
So this looks at how the cost changes as we increase the number of beads
per unit cube, as the number of non-bonded interactions increases with
density.

### oilwater_{angles,hookean,monomer}

These are all a random mixture of water and an oil. For the angles
and hookean versions the oil consists of polymers of length 6, and
for the angle version the polymers are also straight. Apart from the
presences or obsence of bonds, the simulations are the same.

### Bilayer

A Water/Amphiphile bilayer. A bunch of stiff lipids arranges into a
membrance. The height of the simulation is always 32, and the width
and depth are scaled. That means the ratio of lipids to water stays
constant as the volume scales.

## Membrane_plus_polymers

A membrane consisting of two lipids, one common, the other uncommon.
The uncommon one has a head group that is attractive to another set
of free-floating polymers randomly scattered through the volume.
This is one of the three experiments considered in the membranes
journal paper. Volume is scaled by holding height constant and
varying width and depth, so the polymer ratios stay the same.

Files
-----

The files are presented in two forms:

`states_bin` : The raw simulation states, encoded in poets-dpd (binary) format

`xml` : POETS XMLv4 self-contained applications. These have two variants:

    - bonds : Includes full support for hookean and angle-bonds
    - nobonds : hookean and angle bond handlers are turned off. Small performance boost, big reduction in code size.