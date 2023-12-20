This code uses a SPH kernel to convert a particle distribution to a uniform grid representation. F2py is used to create an importable Python module.

INPUT:
- ncores --------------> How many CPUs OMP will use
- partNum -------------> number of particles
- x_pos, y_pos, z_pos -> particle positions
- field ---------------> particle field that is going to be converted to a uniform grid
- Lx, Ly, Lz ----------> box size
- nx, ny, nz ----------> number of cells in each direction for the uniform grid
- kneigh --------------> $k$-nearest neighbour to calculate the $h$ kernel distance.

OUTPUT:
- grid_field ----------> particle field defined on the grid
- hpart ---------------> SPH $h$ kernel distance of each particle

