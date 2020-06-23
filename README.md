# randomparticles

This package is designed to mimic particles distributed in a random particle film for photonics calculations. 

This package generates particles, with radii specified by a set distribution, in random spatial locations over a rectangular region. 
The code can utilize periodic boundary conditions along the X and Y directions to mimic an infinite film, or use hard spherical boundary conditions.
Input parameters are the dimensions of the region, particle radii distribution, and film fill fraction. 

Limitations to the code:
The dimensions of the rectangular region are determined by the nearest integer unit cell to a FCC lattice. 
The code may fail to create fill fractions near the theoretical limits of sphere packing for radii that are not all the same. 





