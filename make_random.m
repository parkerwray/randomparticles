function cords = make_random(cords, r, lower_bounds, upper_bounds, giggles, dimension)

% This function randomizes nanoparticles in a pre-defined box. and applies
% periodic boundary conditions

disp('Randomizing the nanoparticle distribution.')
disp('Itteration number:')
lineLength = fprintf(num2str(1),21.1);
times_moved = zeros(size(cords,1),1);
for i = 1:giggles
    if ~mod(i,10)
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf(num2str(i),21.1);
        pause(0.001);
    end
    for idx = 1:size(cords,1)
       
        if cords(idx,1) == 0 && cords(idx,2) == 0 && cords(idx,3) == 0 
            continue; % dont move the [0,0,0] particle!
        end
        
        new_cord = random_move(cords(idx, :), dimension);
        %new_cord = periodic_BC_3D(new_cord, lower_bounds, upper_bounds);
        new_cord0 = periodic_BC_3D(new_cord, lower_bounds, upper_bounds);
        new_cord = make_mirror_3D(new_cord0, r, lower_bounds, upper_bounds);
        
        dummy = cords;
        dummy(idx,:) = [];
        %flag = check_touch(new_cord, dummy, r, lower_bounds, upper_bounds);
        
        for idx2 = 1:size(new_cord,1)
            flag = check_touch(new_cord(idx2,:), dummy, r, lower_bounds, upper_bounds);
        end
        
        if flag == 0
            times_moved(idx) = times_moved(idx)+1;
            cords(idx,:) = new_cord0;
        end


    end

end
figure,
bar(times_moved)
xlabel('Particle Number')
ylabel('Times Moved')
disp(newline)
end
