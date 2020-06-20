function [new_radii, new_cords] = full_randomize(cords, radii, bounds, ...
    giggles, dimension)

new_cords = cords;
disp(['Randomly removing particles from FCC lattice to match FF request.'])
disp(['Removing ', num2str(size(new_cords, 1)-length(radii)), ' particles.'])


while size(new_cords, 1) > length(radii) % Size function must have dimension argument even for vector.
   % Do not allow removal of [0,0,0] in index = 1;
   new_cords(randi([2,size(new_cords, 1)]),:)=[]; 
end

new_cords = fix_overlap(new_cords, radii, bounds(1,:), bounds(2,:));
new_cords = make_random(new_cords, radii, bounds(1,:), ...
    bounds(2,:), giggles, dimension);
[new_radii, new_cords] = make_all_mirrors(new_cords, radii, ...
    bounds(1,:), bounds(2,:));


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

end