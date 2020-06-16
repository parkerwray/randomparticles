function new_loc = periodic_BC_1D(loc, lower_bound, upper_bound)
% Apply periodic boundary conditions to a 1D vector (used as a constructor
% for 3D vectors
    new_loc = loc;
    len = abs(upper_bound-lower_bound);
    if loc < lower_bound
        new_loc = loc+len;
    elseif loc >= upper_bound
        new_loc = loc-len;
    end

end
