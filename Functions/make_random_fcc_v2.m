function [radii, cord, bounds, a, am] = ...
    make_random_fcc_v2(r, radii, ff, bounds, Nspheres, giggles, dimension)

% NOTE MAX FF IS 55

upper_bound = bounds;
lower_bound = -bounds;
lower_bound(3) = 0;


Vsphere = (4/3)*pi*(r^3);
Asphere = pi*(r^2);
a = 4*r./sqrt(2);  % sqrt(2)*a/2 = 2*r
am = a;

% disp(['Lattice spacing: ', num2str(am)]);
% disp(['Particle spacing : ', num2str(2*d)]);

if upper_bound(3)*r <= r || dimension == 2
    dimension = 2;
    upper_bound(3) = 1.1;
else 
    dimension = 3;
end

if dimension == 2
    % Primative lattice vectors for the NP positions in FCC
    l0 = [0,0,0];
    l1 = [am,am,0]./2;
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
% % Translation vectors to move the lattice across space
% vx = [1,0,0].*am;
% vy = [0,1,0].*am;
% vz = [0,0,1].*am;

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
                %if dimension == 3
                    cord = [cord; l3+v];
                %end
            end
            
        end
    end
end
cord = unique(round(cord,5),'rows');
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

% This variable is for checking degree of randomness vs. periodic.
cord_periodic = cord;

% Determine the number of spheres needed to satisfy fill fraction
% requirement
Vspace = prod((2.*bounds).*a);
Aspace = prod((2.*bounds(1:2)).*a);

%if dimension == 2
%    Nspheres = floor(ff*Aspace/Asphere);
%else
%    Nspheres = floor(ff*Vspace/Vsphere); % Get number of spheres for true ff
%end

%disp(['Requested fill fraction: ', num2str(100*ff)]);
%disp(['Requires ' , num2str(Nspheres), ' particles.']);

Nsize = size(cord,1);
while Nsize > Nspheres
    % Pick a random particle (Not the origin particle) and remove it until
    % you reach Nspheres 
    cord(randi([2,Nsize]),:) = [];
    Nsize = size(cord,1);
end

%removed_periodic_cords = cord;

%Nspheres = size(cord,1);
%if dimension == 2
%    disp(['Created area fill fraciton: ', num2str(100*Nspheres.*Asphere/Aspace)]);
%else
%    disp(['Created fill fraction: ', num2str(100*Nspheres.*Vsphere./Vspace)]);
%end
%disp(['Using ' , num2str(Nspheres), ' particles.']);

bounds = [-upper_bound;upper_bound].*a;
% Add separability part
cord = fix_overlap(cord, radii, bounds(1,:), bounds(2,:));
% Randomize positions after separation
cord = make_random(cord, radii, bounds(1,:), bounds(2,:), giggles, dimension);

% Generates all necessary mirror images for it to be periodic
[radii, cord] = make_all_mirrors(cord, radii, bounds(1,:), bounds(2,:));


disp(newline);
end

