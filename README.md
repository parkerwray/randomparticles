# randomparticles

This package is designed to mimic particles distributed in a random particle film for photonics calculations. 

This package generates particles, with radii specified by a set distribution, in random spatial locations over a rectangular region. 
The code utilizes periodic boundary conditions along the X and Y directions to mimic an infinate film.
Input parameters are the dimensions of the rectangular region, particle radii distribution, and film fill fraction. 

Limitations to the code:
The dimensions of the rectangular region are determined by the nearest interger unit cell to a FCC lattice. 
The code may fail to create fill fractions near the theoretical limits of sphere packing for radii that are not all the same. 





