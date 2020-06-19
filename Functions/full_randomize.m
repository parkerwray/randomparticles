function [new_radii, new_cords] = full_randomize(cords, radii, bounds, ...
    giggles, dimension)

new_cords = cords;
while size(new_cords, 1) > size(radii)
   new_cords(randi(size(new_cords, 1)),:)=[]; 
end
%disp("b");
new_cords = fix_overlap(new_cords, radii, bounds(1,:), bounds(2,:));
%disp("c");
new_cords = make_random(new_cords, radii, bounds(1,:), ...
    bounds(2,:), giggles, dimension);
%disp("d");
[new_radii, new_cords] = make_all_mirrors(new_cords, radii, ...
    bounds(1,:), bounds(2,:));
end