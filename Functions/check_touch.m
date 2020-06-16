function flag = check_touch(new_cord, cords, r, lower_bounds, upper_bounds)

% This function determines if a particle (centered at new_cord, with radius
% r) is touching any other particles (given by particle list cords, and
% radius r)

flag = 0;
for idx = 1:size(cords,1)
    dist = get_distance(new_cord, cords(idx, :), lower_bounds, upper_bounds);
    if dist <= 2*r
        flag = 1;
        break;
    end
end

end