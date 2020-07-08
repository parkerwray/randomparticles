function cords = make_random(cords, r, lower_bounds, upper_bounds, giggles, dimension, loud)

% This function randomizes nanoparticles in a pre-defined box. and applies
% periodic boundary conditions
if nargin < 7
    loud = 1;
end
if loud
    disp('Randomizing the nanoparticle distribution.')
    disp('Iteration number:')
    lineLength = fprintf(num2str(1),21.1);
end
times_moved = zeros(size(cords,1),1);
dist = times_moved;
adist = dist;
original_cords = cords;
for i = 1:giggles
    if loud && ~mod(i,10)
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf(num2str(i),21.1);
        pause(0.001);
    end
    for idx = 1:size(cords,1)
       
        if cords(idx,1) == 0 && cords(idx,2) == 0 && cords(idx,3) == 0 
            continue; % dont move the [0,0,0] particle!
        end
        
        new_cord = random_move(cords(idx, :),r(idx), dimension);
        new_cord0 = periodic_BC_3D(new_cord, lower_bounds + ...
            [0,0,r(idx)], upper_bounds - [0,0,r(idx)]);
        new_cord = make_mirror_3D(new_cord0, r(idx), ...
            lower_bounds, upper_bounds);
        
        dummy = cords;
        dummy(idx,:) = [];
        dummyr = r;
        dummyr(idx) = [];
        %flag = check_touch(new_cord, dummy, r, lower_bounds, upper_bounds);
        flag = 0;
        for idx2 = 1:size(new_cord,1)
            flag = flag + check_touch(new_cord(idx2,:), dummy, r(idx), ...
                dummyr, lower_bounds, upper_bounds);
        end
        %flag = check_touch(new_cord0, dummy, r(idx), dummyr, ...
        %    lower_bounds, upper_bounds);
        if flag == 0
            times_moved(idx) = times_moved(idx)+1;
            [d, ~] = get_distance(cords(idx,:), new_cord0,...
                    lower_bounds, upper_bounds);
            cords(idx,:) = new_cord0;
            dist(idx) = dist(idx)+d;
        end


    end

end

for idx = 1:size(cords,1)
    [adist(idx), ~] = get_distance(cords(idx,:), original_cords(idx,:),...
            lower_bounds, upper_bounds);
end
if loud
figure,
subplot(1,3,1)
bar(times_moved)
xlabel('Particle Number')
ylabel('Times Moved')
subplot(1,3,2)
bar(dist)
xlabel('Particle Number')
ylabel('Total Distance Traveled')
subplot(1,3,3)
bar(adist)
xlabel('Particle Number')
ylabel('Actual Distance Traveled')

disp(newline)
end
end