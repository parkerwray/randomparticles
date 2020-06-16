function [cord, bounds, a, am, Nspheres] = ...
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


% Vsphere = (4/3)*pi*(r^3);
% Asphere = pi*(r^2);
a = 4*r./sqrt(2); 
am = a;

% disp(['Lattice spacing: ', num2str(am)]);
% disp(['Particle spacing : ', num2str(2*d)]);


if dimension == 2
    % Primative lattice vectors for the NP positions in FCC
    l0 = [0,0,0];
    l1 = [am,am,0];
    l2 = [0,am,0];
    l3 = [am,0,0];
    % Translation vectors to move the lattice across space
    vx = [1,0,0].*am;
    vy = [0,1,0].*am;
    vz = [0,0,0].*am;    
else
    % Primative lattice vectors for the NP positions in FCC
    l0 = [0,0,0];
    l1 = [am,am,0]./2;
    l2 = [0,am,am]./2;
    l3 = [am,0,am]./2;
    % Translation vectors to move the lattice across space
    vx = [1,0,0].*am;
    vy = [0,1,0].*am;
    vz = [0,0,1].*am;    
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
  
cord(cord(:,1)<lower_bound(1).*a,:) = [];
cord(cord(:,2)<lower_bound(2).*a,:) = [];
cord(cord(:,3)<lower_bound(3).*a,:) = [];

cord(cord(:,1)>upper_bound(1).*a,:) = [];
cord(cord(:,2)>upper_bound(2).*a,:) = [];
cord(cord(:,3)>upper_bound(3).*a-r,:) = [];

cord= do_z_flip(cord);

% Account for periodic repetitions! 
cord(cord(:,1)>upper_bound(1).*a-r,:) = [];
cord(cord(:,2)>upper_bound(2).*a-r,:) = [];
cord(cord(:,3)>upper_bound(3).*a-r,:) = [];




end