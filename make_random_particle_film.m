


function [spheres, ff] = make_random_particle_film(dimension, boundary_type, r, probability, ff_desired, center_radius, size, giggles, loud, sim_folder, save_folder)

%{
    This function randomly places spheres in a 2D or 3D geomeery where the
    center particle is deterministically placed at the origin and has a
    deterministically (user set) radius.

    The distribution of particle radius is set by a user defines 
    probability distribution determined by the variables, r and
    distribution.
    
    The spatial distribution of particle placement is set by a uniform
    distribution.

    This code has been tested for bugs and has multiple self consistency
    checks implemented at runtime. If an error is caught, the code goes
    into debug mode at that locaiton. 

%} 

addpath(genpath(strcat(sim_folder))); %Add all files in the nanoparticle folder so sim can run.

distr = @(~) randsample(r, 1, true, probability);


margin = 0.01;
tic;
if strcmp(boundary_type, "sphere") == 1
    
    scale = 1;
    
    [radii, ff, Nspheres] = get_radii_and_ff_in_sphere(scale, size, ...
        center_radius, ff_desired, distr, margin, dimension, loud);
    
    % If you are simulating particles based on a spherical boundary then
    % size is the radius of that boundary (i.e., the radius of the largest
    % sphere that will circumscribe all the particles inside). 
    
    cords = full_randomize_in_sphere(radii, size, giggles, dimension, loud);
    radii = radii.';


elseif strcmp(boundary_type, "film") == 1
    
    % If you are simulating a particle film, then the size of the film is a
    % vector [x max, y max, z max]. This is then multiplied by the lattice
    % spaceing abd bounded by both sized. The final dimension would be 
    % 2a*[x max, y max, z max].
    
    
    bounds = size;
    if dimension == 3
        [cords, bounds, a] = make_fcc_3D(min(r), bounds, loud);
    else
        bounds(3) = 1000; % A patch because if you dont make the z-bounds large you will make mirrors of big particles. bec. FCC is dactated by min(r).
        [cords, bounds, a] = make_fcc_2D(min(r), bounds, loud);
        
    end
    [radii, ff, Nspheres] = get_radii_and_ff(bounds, a,...
        center_radius, ff_desired, distr, margin, dimension, loud);
    [radii, cords] = full_randomize_v2(cords, radii, bounds.*a, ...
        giggles, dimension, ff_desired, margin, loud);
else
    disp("Invalid geometry requested."); 
    keyboard;
end

if loud
figure, 
plot_radii(radii); %pass
end


    
%%
if loud
if strcmp(boundary_type, "film") == 1
    make_spheres(cords, radii, bounds(1,:).*a, bounds(2,:).*a);
    disp('Simulation region:')
    disp(['X distance: ', num2str((bounds(2,1)-bounds(1,1)).*a)])
    disp(['Y distance: ', num2str((bounds(2,2)-bounds(1,2)).*a)])
    disp(['Z distance: ', num2str((bounds(2,3)-bounds(1,3)).*a)])
elseif strcmp(boundary_type, "sphere") == 1
    make_spheres_in_sphere(cords, radii, size);
    disp(['Cluster diameter: ',num2str((2*size))])
    %make_spheres(cords, radii, [-1000,-1000,-1000], [1000,1000,1000]);
end
end

spheres = [radii, cords];
spheres(abs(spheres)<(10^-3)) = 0;

if exist('save_folder')
    cd(save_folder)
    writematrix(spheres)    
end
    
end






















