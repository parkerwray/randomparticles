function make_random_particle_distribution(cord, bounds, ff, a, giggles, dimension) 



    make_random_particle_distribution.make_random =...
        @make_random;







    % Determine the number of spheres needed to satisfy fill fraction
    % requirement
    Vspace = prod((2.*bounds).*a);
    Aspace = prod((2.*bounds(1:2)).*a);

    if dimension == 2
        Nspheres = floor(ff*Aspace/Asphere);
    else
        Nspheres = floor(ff*Vspace/Vsphere); % Get number of spheres for true ff
    end

    disp(['Requested fill fraction: ', num2str(100*ff)]);
    disp(['Requires ' , num2str(Nspheres), ' particles.']);

    Nsize = size(cord,1);
    while Nsize > Nspheres
        % Pick a random particle (Not the origin particle) and remove it until
        % you reach Nspheres 
        cord(randi([2,Nsize]),:) = [];
        Nsize = size(cord,1);
    end

    Nspheres = size(cord,1);
    if dimension == 2
        disp(['Created area fill fraciton: ', num2str(100*Nspheres.*Asphere/Aspace)]);
    else
        disp(['Created fill fraction: ', num2str(100*Nspheres.*Vsphere./Vspace)]);
    end
    disp(['Using ' , num2str(Nspheres), ' particles.']);

    bounds = [-upper_bound;upper_bound].*a;

    cord = make_random(cord, r, bounds(1,:)+[0,0,r], bounds(2,:)-[0,0,r], giggles, dimension);

    %cord = make_all_mirrors(cord, r, bounds(1,:), bounds(2,:));

    disp(newline);
end

%% Subfunctions

function cords = make_random(cords, r, lower_bounds, upper_bounds)

% This function randomizes nanoparticles in a pre-defined box. and applies
% periodic boundary conditions

make_random.random_move = @random_move;


disp('Randomizing the nanoparticle distribution.')
disp('Itteration number:')
lineLength = fprintf(num2str(1),21.1);
times_moved = zeros(size(cords,1),1);
for i = 1:2000
    if ~mod(i,10)
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf(num2str(i),21.1);
        pause(0.001);
    end
    for idx = 1:size(cords,1)
       
        if cords(idx,1) == 0 && cords(idx,2) == 0 && cords(idx,3) == 0 
            continue; % dont move the [0,0,0] particle!
        end
        
        new_cord = make_random.random_move(cords(idx, :));
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

function new_cord = random_move(cord)
% Randomly move the particle +1, 0, or -1 in the x, y, and z dimensions.


   new_cord = [cord(1)+round(2*rand-1),...
       cord(2)+round(2*rand-1),...
       cord(3)+round(2*rand-1)];
  
end

function new_cord = periodic_BC_3D(cord, lower_bounds, upper_bounds)
% Apply periodic boundary conditions to a 3D vector.
    xloc = periodic_BC_1D(cord(1), lower_bounds(1), upper_bounds(1));
    yloc = periodic_BC_1D(cord(2), lower_bounds(2), upper_bounds(2));   
    zloc = periodic_BC_1D(cord(3), lower_bounds(3), upper_bounds(3));
    
    new_cord = [xloc, yloc, zloc];
end


