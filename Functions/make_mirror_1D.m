function loc = make_mirror_1D(loc, r, lower_bound, upper_bound)
% This function gives the mirror image coordinate for a particle that
% overlaps with a boundary for periodic BCs

len = abs(upper_bound-lower_bound);

if loc+r > upper_bound
    loc = loc-len;
elseif loc-r < lower_bound
    loc = loc+len;
end

end
