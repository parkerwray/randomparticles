function [cord, bounds, a] = ...
    make_fcc_v3(r, bounds, dimension)
%{ 
    This function generates spheres in a FCC orientation within the limits
    of rectangular bounds. Periodic boundary conditions are used. If a
    sphere lies on an edge of the bounds, it would require a mirror sphere
    to represent the periodic overlap. This mirror sphere is not saved into
    the list of sphere. I.e., only one sphere is saved for spheres at the
    bounds. Generating the mirror sphere can be achieved by a seperate
    function. 
%}


% NOTE MAX FF IS 55

upper_bound = bounds;
lower_bound = -bounds;
lower_bound(3) = 0;

a = 4*r./sqrt(2); 

if dimension == 2
    % Primative lattice vectors for the NP positions in FCC
    l0 = [0,0,0];
    l1 = [a,a,0];
    l2 = [0,a,0];
    l3 = [a,0,0];
    % Translation vectors to move the lattice across space
    vx = [1,0,0].*a;
    vy = [0,1,0].*a;
    vz = [0,0,0].*a;    
else
    % Primative lattice vectors for the NP positions in FCC
    l0 = [0,0,0];
    l1 = [a,a,0]./2;
    l2 = [0,a,a]./2;
    l3 = [a,0,a]./2;
    % Translation vectors to move the lattice across space
    vx = [1,0,0].*a;
    vy = [0,1,0].*a;
    vz = [0,0,1].*a;    
end


% Generate 
cord = [0,0,0];
v = [0,0,0];
for z = lower_bound(3):upper_bound(3)   
    for y = lower_bound(2):upper_bound(2)
        vzz = z.*vz;
        vyy = vzz + y.*vy;
        
        for x = lower_bound(1):upper_bound(1)
            v = vyy+x.*vx;
            
            if x==0 && y==0 && z==0
                continue;
            else
                cord = [cord; l0+v];
                cord = [cord; l1+v];
                cord = [cord; l2+v];
                cord = [cord; l3+v];
            end
            
        end
    end
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

end