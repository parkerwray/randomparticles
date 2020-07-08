function new_cords = full_randomize_in_sphere(radii, ...
    R, giggles, dimension, loud)

new_cords = zeros(length(radii), 3);

new_cords(1,:) = [0,0,0];

for idx = 2:size(new_cords,1)
    r = rand() * (R - radii(idx));
    theta = rand() * 2 * pi;
    if dimension == 3
        phi = rand() * pi;
    else
        phi = pi / 2;
    end
    new_cords(idx, :) = [r*cos(theta)*sin(phi), r*sin(theta)*sin(phi),...
        r*cos(phi)];
end


new_cords = fix_overlap_sphere(new_cords, radii, R, loud);


new_cords = make_random_sphere(new_cords, radii, R, giggles, dimension, loud);

%{
FLAG = check_zero_repeat_overlap(new_cords, bounds, new_radii);
    
if FLAG == 0
    disp('Full randomization sucessful!')
    disp('A particle is at [0,0,0].')
    disp('No repeat particles were found.')
    disp('No overlaping particles were found.')
    disp(newline)
else
    keyboard;
end
%}

end