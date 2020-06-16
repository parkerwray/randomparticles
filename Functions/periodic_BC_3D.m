function new_cord = periodic_BC_3D(cord, lower_bounds, upper_bounds)
% Apply periodic boundary conditions to a 3D vector.
    xloc = periodic_BC_1D(cord(1), lower_bounds(1), upper_bounds(1));
    yloc = periodic_BC_1D(cord(2), lower_bounds(2), upper_bounds(2));   
    zloc = periodic_BC_1D(cord(3), lower_bounds(3), upper_bounds(3));
    
    new_cord = [xloc, yloc, zloc];
end

