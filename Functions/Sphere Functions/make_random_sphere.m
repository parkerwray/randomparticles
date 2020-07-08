function cords = make_random_sphere(cords, r, R, giggles, dimension, loud)

% This function randomizes nanoparticles in a pre-defined sphere. and applies
% periodic boundary conditions
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
    if loud
        if ~mod(i,10)
            fprintf(repmat('\b',1,lineLength))
            lineLength = fprintf(num2str(i),21.1);
            pause(0.001);
        end
    end
    for idx = 1:size(cords,1)
       
        if cords(idx,1) == 0 && cords(idx,2) == 0 && cords(idx,3) == 0 
            continue; % dont move the [0,0,0] particle!
        end
        
        new_cord0 = random_move(cords(idx,:),r(idx), dimension);
        %{
        dummy = cords;
        dummy(idx,:) = [];
        dummyr = r;
        dummyr(idx) = [];
        %}
        flag = 0;
        for idx2 = 1:size(cords, 1)
           if idx ~= idx2
               if norm(new_cord0-cords(idx2,:)) < r(idx)+r(idx2)
                  flag = 1;
                  break;
               end
           end
        end
        if (flag == 0) && (norm(new_cord0) + r(idx) <= R)
            times_moved(idx) = times_moved(idx)+1;
            d = norm(new_cord0 - cords(idx,:));
            cords(idx,:) = new_cord0;
            dist(idx) = dist(idx)+d;
        end


    end

end

for idx = 1:size(cords,1)
    [adist(idx), ~] = get_distance(cords(idx,:), original_cords(idx,:),...
            [-Inf, -Inf, -Inf], [Inf, Inf, Inf]);
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