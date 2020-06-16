
function new_cord = make_mirror_3D(cord, r, lower_bounds, upper_bounds)
% This function gives all the coordinates of a nanoparticle in a square
% periodic BC condition. I.e., it determines the particles center location
% and generates the necessary "mirror image" particles to apply periodic
% BCs to particles on the edges. 

    x0 = cord(1);
    y0 = cord(2);
    z0 = cord(3);

% Apply periodic boundary conditions to a 3D vector.
    xloc = make_mirror_1D(cord(1), r, lower_bounds(1), upper_bounds(1));
    yloc = make_mirror_1D(cord(2), r, lower_bounds(2), upper_bounds(2));   
    zloc = make_mirror_1D(cord(3), r, lower_bounds(3), upper_bounds(3));
    
    new_cord = [x0, y0, z0;...
                xloc, y0, z0;...
                xloc, yloc, z0;...
                xloc, yloc, zloc;...
                x0, yloc, z0;...
                x0, yloc, zloc;...
                x0, y0, zloc;...
                xloc, y0, zloc];
                
     new_cord = unique(new_cord,'rows');    

end
