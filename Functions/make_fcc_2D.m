function [cord, bounds, a, FLAG] = ...
    make_fcc_2D(r, bounds)
%{ 
    This function generates spheres in a FCC orientation within the limits
    of rectangular bounds. Periodic boundary conditions are used. If a
    sphere lies on an edge of the bounds, it would require a mirror sphere
    to represent the periodic overlap. This mirror sphere is not saved into
    the list of sphere. I.e., only one sphere is saved for spheres at the
    bounds. Generating the mirror sphere can be achieved by a seperate
    function. 
%}
FLAG = 0;
upper_bound = bounds;
lower_bound = -bounds;
lower_bound(3) = 0;

a = 4*r./sqrt(2); 

% Primative lattice vectors for the NP positions in FCC
l0 = [0,0,0];
l1 = [a,a,0]./2;

% Translation vectors to move the lattice across space
vx = [1,0,0].*a;
vy = [0,1,0].*a;  

% Generate particles in 2D FCC 
%cord = [0,0,0];
cord = [];
v = [0,0,0];
for y = lower_bound(2):upper_bound(2)
    vyy = y.*vy;
    for x = lower_bound(1):upper_bound(1)
        v = vyy+x.*vx;
        cord = [cord; l0+v];
        cord = [cord; l1+v];
        % Remove comment to visualize for debug.
%         r2 = r.*ones(size(cord,1),1); %account for all same radius using updated functions
%         bounds = [-upper_bound; upper_bound];
%         make_spheres_in_rectangle(cord, r2, bounds(1,:).*a, bounds(2,:).*a);
    end  
end

% Validate that there are no particle repeats (i.e., complete overlap).
% If there are particle repeats, remove all but one of the repeat particles
% and send an error to the user that the simulation may be incorrect. 
dummy = unique(cord,'rows'); %Breaks 3D necessary for 2D
if size(dummy,1) ~= size(cord,1)
    FLAG = 1;
    disp('ERROR! When making the FCC lattice particles were repeated.')
    disp('Repeat particles are being removed.')
    disp('Fill fraction and other metrics may be incorrect!')
    cord = dummy;
end
[flag,loc] = ismember([0,0,0],cord,'rows');
if flag == 1
    cord(loc,:) = [];
    cord = [0,0,0;cord];
else
    FLAG = 1;
    disp('ERROR! No particle was generated in the center!')
    %keyboard;
end

% An FCC lattice is generated slightly overflowing the desired simulation
% region. This is becuase you repeat an enire lattice, which has multiple
% particles. Therefore, to get the simulation region we want, we cut our
% the particles that are over our simulation bounds. 

% Cut particles in lower half of x,y,z bounds  
cord(cord(:,1)<lower_bound(1).*a,:) = [];
cord(cord(:,2)<lower_bound(2).*a,:) = [];
cord(cord(:,3)<lower_bound(3).*a,:) = [];

% Cut particles in lower half of x,y,z bounds 
cord(cord(:,1)>upper_bound(1).*a,:) = [];
cord(cord(:,2)>upper_bound(2).*a,:) = [];
cord(cord(:,3)>upper_bound(3).*a-r,:) = [];

% To ensure that the center particle is always at the center of the entire 
% simulation region, the FCC lattice is generated first only in the upper
% half plane z = 0 to z_max. Then the entire particle distribution is 
% flipped onto the lower half plane.
cord= do_z_flip(cord);

% Account for periodic repetitions! Only keep one particle for every
% particle on the edge of the simulation region.
cord(cord(:,1)>upper_bound(1).*a-r,:) = [];
cord(cord(:,2)>upper_bound(2).*a-r,:) = [];
cord(cord(:,3)>upper_bound(3).*a-r,:) = [];

bounds = [-upper_bound; upper_bound];
bounds(:,3) = bounds(:,3).*r/a; % Z length = 2r bec. 2D.

% Check and make sure no particles are overlapping. 
[overlap_distance, overlap_idx] = ...
    check_distance_function_v2(cord, r, bounds(1,:).*a, bounds(2,:).*a);

if sum(overlap_idx) ~= 0
    disp(['ERROR! Particles created in the 2D FCC method are overlaping!'])
    disp(['This may break other functions and return wrong FF metrics'])
    FLAG = 1;
end


if FLAG == 0
    disp('2D FCC lattice generation was sucessful!')
    disp('A particle was generated at [0,0,0].')
    disp('No repeat particles were found.')
    disp('No overlaping particles were found.')
else
    keyboard;
end

end