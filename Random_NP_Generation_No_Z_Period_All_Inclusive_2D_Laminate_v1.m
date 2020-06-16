%% RUN TESTS
clear;
clc;
% test_PBCs; % Pass
% test_fcc; % Pass
% test_flips; % Pass
% test_make_random_v2; % Pass

r = 44;
ff = 0.20;

bounds = [10.5,10.5,0.5]; 
%bounds = [3,3,0.5]; 
giggles = 20000;
dimension  = 2;
tic;

[cords, bounds, a, am, Nspheres] = ...
    make_random_fcc_v2(r, ff, bounds, giggles, dimension);
%%
% clc;
% lower_bound = bounds(1,:);
% upper_bound = bounds(2,:);
% Itop = find_top_view_ff(cords, r, lower_bound, upper_bound);
% fftop = sum(sum(Itop,1),2)/(size(Itop,1).*size(Itop,2));
% figure, imagesc(Itop)
% [Ntouch, ffcalc, is_zero, norm_dist, overlap_idx] = ...
%     check_touch_ff_orig(cords, ff, r, bounds(1,:), bounds(2,:));


%%
make_spheres(cords, r, bounds(1,:), bounds(2,:));
disp(['X distance: ', num2str(bounds(2,1)-bounds(1,1)), '. Multiples of r: ', num2str((bounds(2,1)-bounds(1,1))/r)])
disp(['Y distance: ', num2str(bounds(2,2)-bounds(1,2)), '. Multiples of r: ', num2str((bounds(2,2)-bounds(1,2))/r)])
disp(['Z distance: ', num2str(bounds(2,3)-bounds(1,3)), '. Multiples of r: ', num2str((bounds(2,3)-bounds(1,3))/r)])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTIONS FOR PLOTTING SPHERES IN FANCY WAY
function I = find_top_view_ff(cords, r, lower_bound, upper_bound)
L = upper_bound-lower_bound;
I = zeros(round(L(1)), round(L(2)));


for idx = 1:size(cords,1)
    cord = [cords(idx,1), cords(idx,2)];
    Idummy = make_circle(I, cord, r, lower_bound, upper_bound);
    I = I|Idummy;
end


% for x = 1:size(circle,1)
%     for y = 1:size(circle,2)
%         if sqrt((center-x).^2+(center-y).^2) <= r
%             circle(x,y) = 1;
%         end
%     end
% end
% 
% for idx = 1:size(cords,1)
%     xloc = round(cord(idx,1)+lower_bound(1));
%     yloc = round(cord(idx,2)+lower_bound(2));
% end
    

 
end

function I = make_circle(I, cord, r, lower_bound, upper_bound)

center = [size(I,1)/2,size(I,2)/2];
for idxx = 1:size(I,1)
    for idxy = 1:size(I,2)
        x = idxx-center(1);
        y = idxy-center(2);
        if sqrt((cord(1)-x).^2+(cord(2)-y).^2)<= r
            I(idxx,idxy) = 1;
        end
    end
end

end

function make_spheres(cord, r, lower_bound, upper_bound)

figure, 
transparency = 0.2;
plotcube(upper_bound-lower_bound, lower_bound, transparency);
hold on 

[x,y,z] = sphere;
% figure,
hold on 
for i = 1:size(cord,1)
    surf(r.*x+cord(i,1), r.*y+cord(i,2), r.*z+cord(i,3))
    %shading interp
end
title(['Number of particles: ', num2str(size(cord,1))])

end

function plotcube(varargin)
% PLOTCUBE - Display a 3D-cube in the current axes
%
%   PLOTCUBE(EDGES,ORIGIN,ALPHA,COLOR) displays a 3D-cube in the current axes
%   with the following properties:
%   * EDGES : 3-elements vector that defines the length of cube edges
%   * ORIGIN: 3-elements vector that defines the start point of the cube
%   * ALPHA : scalar that defines the transparency of the cube faces (from 0
%             to 1)
%   * COLOR : 3-elements vector that defines the faces color of the cube
%
% Example:
%   >> plotcube([5 5 5],[ 2  2  2],.8,[1 0 0]);
%   >> plotcube([5 5 5],[10 10 10],.8,[0 1 0]);
%   >> plotcube([5 5 5],[20 20 20],.8,[0 0 1]);
% Default input arguments
inArgs = { ...
  [10 56 100] , ... % Default edge sizes (x,y and z)
  [10 10  10] , ... % Default coordinates of the origin point of the cube
  .7          , ... % Default alpha value for the cube's faces
  [1 0 0]       ... % Default Color for the cube
  };
% Replace default input arguments by input values
inArgs(1:nargin) = varargin;
% Create all variables
[edges,origin,alpha,clr] = deal(inArgs{:});
XYZ = { ...
  [0 0 0 0]  [0 0 1 1]  [0 1 1 0] ; ...
  [1 1 1 1]  [0 0 1 1]  [0 1 1 0] ; ...
  [0 1 1 0]  [0 0 0 0]  [0 0 1 1] ; ...
  [0 1 1 0]  [1 1 1 1]  [0 0 1 1] ; ...
  [0 1 1 0]  [0 0 1 1]  [0 0 0 0] ; ...
  [0 1 1 0]  [0 0 1 1]  [1 1 1 1]   ...
  };
XYZ = mat2cell(...
  cellfun( @(x,y,z) x*y+z , ...
    XYZ , ...
    repmat(mat2cell(edges,1,[1 1 1]),6,1) , ...
    repmat(mat2cell(origin,1,[1 1 1]),6,1) , ...
    'UniformOutput',false), ...
  6,[1 1 1]);
cellfun(@patch,XYZ{1},XYZ{2},XYZ{3},...
  repmat({clr},6,1),...
  repmat({'FaceAlpha'},6,1),...
  repmat({alpha},6,1)...
  );
view(3);


end

%% MAKE FCC PERIODIC ARRAY %%

function [cord, bounds, a, am, Nspheres] = ...
    make_random_fcc_v2(r, ff, bounds, giggles, dimension)

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

cord = make_all_mirrors(cord, r, bounds(1,:), bounds(2,:));


disp(newline);
end


%% MAKE 3D VERSIONS OF PERIODIC BCS %%

function final_cords = make_all_mirrors(cords, r, lower_bound, upper_bound)
    
    final_cords = [];    
    for idx2 = 1:size(cords,1)
        dummy = make_mirror_3D(cords(idx2,:), r, lower_bound, upper_bound);
        final_cords = [final_cords;dummy];
    end


end

function new_cord = make_mirror_3D(cord, r, lower_bounds, upper_bounds)
% This function gives all the coordinates of a nanoparticle in a square
% periodic BC condition. I.e., it determines the particles center location
% and generates the necessary "mirror image" particles to apply periodic
% BCs to particles on the edges. 

    x0 = cord(1);
    y0 = cord(2);
    z0 = cord(3);

% Apply periodic boundary conditions to a 3D vector.
    xloc = make_mirror_1D(cord(1), r, lower_bounds(1), upper_bounds(1));
    yloc = make_mirror_1D(cord(2), r, lower_bounds(2), upper_bounds(2));   
    zloc = make_mirror_1D(cord(3), r, lower_bounds(3), upper_bounds(3));
    
    new_cord = [x0, y0, z0;...
                xloc, y0, z0;...
                xloc, yloc, z0;...
                xloc, yloc, zloc;...
                x0, yloc, z0;...
                x0, yloc, zloc;...
                x0, y0, zloc;...
                xloc, y0, zloc];
                
     new_cord = unique(new_cord,'rows');    

end

function new_cord = periodic_BC_3D(cord, lower_bounds, upper_bounds)
% Apply periodic boundary conditions to a 3D vector.
    xloc = periodic_BC_1D(cord(1), lower_bounds(1), upper_bounds(1));
    yloc = periodic_BC_1D(cord(2), lower_bounds(2), upper_bounds(2));   
    zloc = periodic_BC_1D(cord(3), lower_bounds(3), upper_bounds(3));
    
    new_cord = [xloc, yloc, zloc];
end


%% MAKE 1D VERSIONS OF PERIODIC BCS %%  

function loc = make_mirror_1D(loc, r, lower_bound, upper_bound)
% This function gives the mirror image coordinate for a particle that
% overlaps with a boundary for periodic BCs

len = abs(upper_bound-lower_bound);

if loc+r > upper_bound
    loc = loc-len;
elseif loc-r < lower_bound
    loc = loc+len;
end

end

function new_loc = periodic_BC_1D(loc, lower_bound, upper_bound)
% Apply periodic boundary conditions to a 1D vector (used as a constructor
% for 3D vectors
    new_loc = loc;
    len = abs(upper_bound-lower_bound);
    if loc < lower_bound
        new_loc = loc+len;
    elseif loc >= upper_bound
        new_loc = loc-len;
    end

end


%% FUNCTIONS FOR RANDOMIZING PARTICLES %%

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

function [dist, dcord] = get_distance(cord1, cord2, lower_bounds, upper_bounds)

dcord1 = periodic_BC_3D((cord1-cord2), lower_bounds, upper_bounds);
dist1 = sqrt(dcord1(1)^2+dcord1(2)^2+dcord1(3)^2);

dcord2 = periodic_BC_3D((cord2-cord1), lower_bounds, upper_bounds);
dist2 = sqrt(dcord2(1)^2+dcord2(2)^2+dcord2(3)^2);

if dist2 < dist1
    dist = dist2;
    dcord = dcord2;
else
    dist = dist1;
    dcord = dcord1;
end


end

function flag = check_touch(new_cord, cords, r, lower_bounds, upper_bounds)

% This function determines if a particle (centered at new_cord, with radius
% r) is touching any other particles (given by particle list cords, and
% radius r)

flag = 0;
for idx = 1:size(cords,1)
    dist = get_distance(new_cord, cords(idx, :), lower_bounds, upper_bounds);
    if dist <= 2*r
        flag = 1;
        break;
    end
end

end

function new_cord = random_move(cord, dimension)
% Randomly move the particle +1, 0, or -1 in the x, y, and z dimensions.
if dimension == 2
       new_cord = [cord(1)+round(2*rand-1),...
       cord(2)+round(2*rand-1),...
       cord(3)];
else
   new_cord = [cord(1)+round(2*rand-1),...
       cord(2)+round(2*rand-1),...
       cord(3)+round(2*rand-1)];
end
  
end

function [cords, lower_bounds, upper_bounds] = center_cords(cords, r, lower_bounds, upper_bounds)

dummy = make_all_mirrors(cords, r, lower_bounds, upper_bounds);
lower_bounds = lower_bounds-mean(dummy);
upper_bounds = upper_bounds-mean(dummy);
cords = cords-mean(dummy);

[~, test_idx] = min(sum(abs(cords),2));
cords(test_idx,:) = round(cords(test_idx,:));

dummyx = cords(:,1);
dummyx(abs(dummyx)<10^(-12)) = 0;

dummyy = cords(:,2);
dummyy(abs(dummyy)<10^(-12)) = 0;

dummyz = cords(:,3);
dummyz(abs(dummyz)<10^(-12)) = 0;

cords = [dummyx, dummyy, dummyz];
end

function cords = do_z_flip(cords)


top_cords = cords;
bottom_cords = [1,1,-1].*cords;

cords = [top_cords; bottom_cords];
idx_origin_particle = find(cords(:,1)==0& cords(:,2) == 0 & cords(:,3) == 0);
cords(idx_origin_particle,:) = [];
cords = unique(cords,'rows');   
cords = [[0,0,0];cords];


end

function [cords, lower_bound, upper_bound] = do_flips(cords,  lower_bound, upper_bound)

v1 = [1,0,0]; %Y-Z plane orth unit vector
v2 = [0,1,0]; %X-Z plane orth unit vector
v3 = [0,0,1]; %X-Y plane orth unit vector

for idx = 1:size(cords, 1)
    cords = [cords; cords(idx,:)-2.*dot(v1,cords(idx,:)).*v1];
    cords = [cords; cords(idx,:)-2.*dot(v2,cords(idx,:)).*v2];
    cords = [cords; cords(idx,:)-2.*dot(v1,cords(idx,:)).*v1-2.*dot(v2,cords(idx,:)).*v2];
    cords = [cords; cords(idx,:)-2.*dot(v3,cords(idx,:)).*v3];
    cords = [cords; cords(idx,:)-2.*dot(v3,cords(idx,:)).*v3-2.*dot(v1,cords(idx,:)).*v1];
    cords = [cords; cords(idx,:)-2.*dot(v3,cords(idx,:)).*v3-2.*dot(v2,cords(idx,:)).*v2];   
    cords = [cords; cords(idx,:)-2.*dot(v3,cords(idx,:)).*v3-2.*dot(v1,cords(idx,:)).*v1-2.*dot(v2,cords(idx,:)).*v2];
end
cords = unique(cords,'rows');    

xmin = min(cords(:,1));
ymin = min(cords(:,2));
zmin = min(cords(:,3));

xmax = max(cords(:,1));
ymax = max(cords(:,2));
zmax = max(cords(:,3));

lower_bound = [xmin, ymin, zmin];
upper_bound = [xmax, ymax, zmax];

idx_remove = [];
for idx = 1:size(cords,1)
    if cords(idx,1) == xmax || cords(idx,2) == ymax || cords(idx,3) == zmax
        idx_remove = [idx_remove, idx];
    end
end
    
cords(idx_remove,:) = [];

end

%% FUNCTIONS FOR MSTM %%

function repeat = check_repeat(sphere_cords)

    last = length(sphere_cords);
    repeat = 0;
    for idx=1:last-1
        if isequal(sphere_cords{idx}, sphere_cords{last})
            repeat = 1;
            break;
        end
    end
    
end


%% FUNCTIONS FOR TESTING THINGS %%

function [dist2, overlap_idx] = check_distance_function(cords, r, lower_bounds, upper_bounds)
overlap_idx = [];
for idx = 1:size(cords,1)
    for idx2 = 1:size(cords,1)
        dist(idx, idx2) = get_distance(cords(idx2, :), cords(idx, :), lower_bounds, upper_bounds);
        dist2(idx,idx2) = dist(idx,idx2)./(2*r);
        if dist2(idx,idx2) < 1 && dist2(idx, idx2) > 0
            overlap_idx = [overlap_idx;[idx,idx2]];
        end
    end
end



end


function [overlap_distance, overlap_idx] = check_distance_function_v2(cords, r, lower_bounds, upper_bounds)
overlap_idx = [];
overlap_distance = [];
for idx = 1:size(cords,1)
    for idx2 = 1:size(cords,1)
        dist(idx, idx2) = get_distance(cords(idx2, :), cords(idx, :), lower_bounds, upper_bounds);
        dist2(idx,idx2) = dist(idx,idx2)-(2*r);
        if dist2(idx,idx2) < -1*10^(-10) && dist(idx, idx2) > 10^(-10)
            overlap_idx = [overlap_idx;[idx,idx2]];
            overlap_distance = [overlap_distance; dist2(idx,idx2)];
        end
    end
end



end






function [Ntouch, ffcalc, is_zero, dist, overlap_idx] = check_touch_ff_orig(cord, ff, r, lower_bound, upper_bound)

    dummy = cord;
    [overlap_distance, overlap_idx] = check_distance_function_v2(dummy, r, lower_bound, upper_bound);
    Ntouch = size(overlap_idx,1);
    [I, ffcalc] = ...
        visualize_spheres(dummy,...
        round(upper_bound(1)-lower_bound(1)), round(upper_bound(3)-lower_bound(3)),...
        r, ff, 0);
    figure, 
    imshow(squeeze(I(:,:,5)))
    
    
    dist = sqrt(dummy(:,1).^2+dummy(:,2).^2+dummy(:,3).^2);
    
    is_zero = min(abs(dist));

    disp(['Desired fill fraction: ', num2str(100*ff)])
    disp(['Calculated fill fraction: ', num2str(100*ffcalc)])
    disp(['Number of touching spheres: ', num2str(Ntouch)])
    
    if Ntouch > 0
        disp(['Worst fractional overlap: ', num2str(100.*min(overlap_distance)), '%'])
        disp(['Overlaped spheres: '])
            disp(num2str(overlap_idx))
            keyboard;
    end
    
    disp(['Orig particle distance to orig: ', num2str(is_zero)])
    disp(newline)
end

function test_flips
clear;
clc;
close all;

ff = 0.1:0.01:0.7;

parfor idx = 1:length(ff)
    lower_bound = [0,0,0];
    upper_bound = [2,2,2];
    r = 10;

    [cord, bounds, ~, ~] = make_fcc_v2(r, ff(idx), lower_bound, upper_bound);
    [cord_flip, lower_bound, upper_bound] = do_flips(cord, bounds(1,:), bounds(2,:));
    cord_final = make_all_mirrors(cord_flip, r, lower_bound, upper_bound);

    [Ntouch(idx), ffcalc(idx), is_zero(idx)] = ...
        check_touch_ff_orig(cord_final, ff(idx), r, lower_bound, upper_bound);

    
end

figure, 
plot(ff, ffcalc)
xlabel('FF')
ylabel('FF_{SIM}')


figure, 
plot(ff, Ntouch)
xlabel('FF')
ylabel('Ntouch')

figure, 
plot(ff, is_zero)
xlabel('FF')
ylabel('Zero Cord')

[x,y,z] = sphere(r);
figure,
hold on 
for i = 1:size(cord,1)
    surf(r.*x+cord(i,1), r.*y+cord(i,2), r.*z+cord(i,3))
    shading interp
end
title(['Number of particles: ', num2str(size(cord,1))])

% figure,
% hold on 
% for i = 1:size(cord_flip,1)
%     surf(r.*x+cord_flip(i,1), r.*y+cord_flip(i,2), r.*z+cord_flip(i,3))
%     shading interp
% end
% title(['Number of particles: ', num2str(size(cord_flip,1))])
% 
% figure,
% hold on 
% for i = 1:size(cord_final,1)
%     surf(r.*x+cord_final(i,1), r.*y+cord_final(i,2), r.*z+cord_final(i,3))
%     shading interp
% end
% title(['Number of particles: ', num2str(size(cord_final,1))])




end

function [ff_calc, Ntouch, ff, is_zero] = test_fcc

clear;
clc;
ff=0.1;



ff = 0.1:0.01:0.73;

parfor (idx = 1:length(ff))
    lower_bound = [0,0,0];
    upper_bound = [4,4,4];
    r = 10;
    
    [cords, bounds, ~, ~] = make_fcc_v2(r, ff(idx), lower_bound, upper_bound);
    [cords, lower_bound, upper_bound] = center_cords(cords, r, bounds(1,:), bounds(2,:));
    cord_final = make_all_mirrors(cords, r, lower_bound, upper_bound);

    dummy = cord_final;
    norm_dist = check_distance_function(dummy, r, lower_bound, upper_bound);
    values = norm_dist(norm_dist<1);
    values(values == 0) = [];
    Ntouch(idx) = sum(values<1);
    [~, ff_calc(idx)] = ...
        visualize_spheres(dummy,...
        round(upper_bound(1)), round(upper_bound(3)),...
        r, ff(idx), 0);
    dist = sqrt(dummy(:,1).^2+dummy(:,2).^2+dummy(:,3).^2);
    
    is_zero(idx,:) = min(abs(dist));
    
    
    
end

figure, 
plot(ff, ff_calc)
xlabel('FF')
ylabel('FF_{SIM}')


figure, 
plot(ff, Ntouch)
xlabel('FF')
ylabel('Ntouch')

figure, 
plot(ff, is_zero)
xlabel('FF')
ylabel('Zero Cord')

end

function [ff_calc, Ntouch, ff, is_zero] = test_make_random
close all;
clear;
clc;
ff = 0.3;
n = 20;
%delete(gcp('nocreate'));
%parpool(n)
    for (idx = 1:length(ff))
        lower_bound = [0,0,0];
        upper_bound = [4,4,4];
        r = 10;
        [cords, bounds] = make_fcc_v2(r, ff(idx), lower_bound, upper_bound);
        [cords, lower_bound, upper_bound] = center_cords(cords,r, bounds(1,:), bounds(2,:));
        cords = make_random(cords, r, lower_bound, upper_bound);
        final_cords = make_all_mirrors(cords, r, lower_bound, upper_bound);
        
        
        norm_dist = check_distance_function(final_cords, r, lower_bound, upper_bound);
        values = norm_dist(norm_dist<1);
        values(values == 0) = [];
        Ntouch(idx) = sum(values<1);
        [~, ff_calc(idx)] = ...
            visualize_spheres(final_cords,...
            round(upper_bound(1)), round(upper_bound(3)),...
            r, ff(idx), 0);
        is_zero(idx,:) = min(abs(final_cords));
    end


figure, 
plot(ff, ff_calc)
xlabel('FF')
ylabel('FF_{SIM}')


figure, 
plot(ff, Ntouch)
xlabel('FF')
ylabel('Ntouch')

figure, 
plot(ff, is_zero)
xlabel('FF')
ylabel('Zero Cord')


end

function test_make_random_v2

clear;
clc;
close all;

ff = 0.1:0.01:0.7;

parfor idx = 1:length(ff)
    lower_bound = [0,0,0];
    upper_bound = [3,3,3];
    r = 10;

    [cord, bounds, ~, ~] = make_fcc_v2(r, ff(idx), lower_bound, upper_bound);
    [cord_flip, lower_bound, upper_bound] = do_flips(cord, bounds(1,:), bounds(2,:));
    
    cord_final_per = make_all_mirrors(cord_flip, r, lower_bound, upper_bound);
    
    cord_random = make_random(cord_flip, r, lower_bound, upper_bound);
    cord_final = make_all_mirrors(cord_random, r, lower_bound, upper_bound);

    [Ntouch(idx), ffcalc(idx), is_zero(idx)] = ...
        check_touch_ff_orig(cord_final_per, ff(idx), r, lower_bound, upper_bound);

    [Ntouch(idx), ffcalc(idx), is_zero(idx)] = ...
        check_touch_ff_orig(cord_final, ff(idx), r, lower_bound, upper_bound);

    
end

figure, 
plot(ff, ffcalc)
xlabel('FF')
ylabel('FF_{SIM}')


figure, 
plot(ff, Ntouch)
xlabel('FF')
ylabel('Ntouch')

figure, 
plot(ff, is_zero)
xlabel('FF')
ylabel('Zero Cord')

% [x,y,z] = sphere(r);


% figure,
% hold on 
% for i = 25:26 %size(cord_random,1)
%     surf(r.*x+cord_random(i,1), r.*y+cord_random(i,2), r.*z+cord_random(i,3))
%     shading interp
% end
% title(['Number of particles: ', num2str(size(cord_random,1))])

% figure,
% hold on 
% for i = 1:size(cord_final,1)
%     surf(r.*x+cord_final(i,1), r.*y+cord_final(i,2), r.*z+cord_final(i,3))
%     shading interp
% end
% title(['Number of particles: ', num2str(size(cord_final,1))])

% for idx = 1:size(overlap_idx,1)
%     disp(['Particle 1: ', num2str(cord_final(overlap_idx(idx,1),:))])
%     disp(['Particle 2: ', num2str(cord_final(overlap_idx(idx,2),:))])
%     disp(num2str(norm_dist(overlap_idx(idx,1),overlap_idx(idx,2))))
% 
% end


end

function test_PBCs


r = 20;
[x,y,z] = sphere(r);
x = r.*x;
y = r.*y;
z = r.*z;

lower_bounds = [-2,-2,-2].*r;
upper_bounds = -lower_bounds;

sim_size = upper_bounds-lower_bounds;

cord = [sim_size(1)/2, sim_size(2)/2, sim_size(3)/2];
for i = 1:200
    
    cord = cord+[1,1,1];
    cord = periodic_BC_3D(cord, lower_bounds, upper_bounds);
    new_cord = make_mirror_3D(cord, r, lower_bounds, upper_bounds);
    
    clf;
    for j = 1:size(new_cord,1)  
        hold on 

        surf(x+new_cord(j,1), y+new_cord(j,2), z+new_cord(j,3))
        shading interp
        xlim([lower_bounds(1),upper_bounds(1)]);
        ylim([lower_bounds(2),upper_bounds(2)]);
        zlim([lower_bounds(2),upper_bounds(3)]);
        view([-10,-40,30])

    end

    F(i) = getframe;

end
movie(F,500,10)



end



%% JUNK FUNCTIONS %%

function [I, ffcalc] = visualize_spheres(cords, sim_length, sim_height, r, ff, plots)

    I = zeros(sim_length,sim_length, sim_height);
    cords = cords-[min(cords(:,1)),min(cords(:,2)), min(cords(:,3))];
    I = draw_spheres(I, cords, r);
    I(I>1) = 1;

    disp(['Requested fill fraction ', num2str(100*ff), ' %'])
    disp(['Simulated fill fraction ',...
        num2str(100.*sum(I(:))./(size(I,1).*size(I,2).*size(I,3))), '%'])
    ffcalc = sum(I(:))./(size(I,1).*size(I,2).*size(I,3));
    
%     x_center = min(abs(cords(:,1)));
%     y_center = min(abs(cords(:,2)));
%     z_center = min(abs(cords(:,3)));
%     
%     disp(['Center Particle is x = ',num2str(x_center),...
%         ', y = ',num2str(y_center),...
%         ', z = ', num2str(z_center)]);
%        
if plots == 1    
    figure, 
    isosurface(I,0)
    %imshow(I)
end

end

function I2 = draw_spheres(I, cords, r)
dummy = I;
imageSizeX = size(I,1);
imageSizeY = size(I,2);
imageSizeZ = size(I,3);

[XInImage, YInImage, ZInImage] = meshgrid(1:imageSizeX, 1:imageSizeY, 1:imageSizeZ);
% Next create the circle in the image.
radius = r;

for idx1 = 1:size(cords, 1)
    
    circle = [];
    centerX = cords(idx1, 1);
    centerY = cords(idx1, 2);
    centerZ = cords(idx1, 3);
    
    circle= (YInImage - centerY).^2 ...
    + (XInImage - centerX).^2 + (ZInImage - centerZ).^2 <= radius.^2;
        
    dummy = dummy+circle;
    
end
I2 = dummy;
end

function circle = make_sphere(I, cord, r)

imageSizeX = size(I,1);
imageSizeY = size(I,2);
imageSizeZ = size(I,3);

[XInImage, YInImage, ZInImage] = meshgrid(1:imageSizeX, 1:imageSizeY, 1:imageSizeZ);
% Next create the circle in the image.
centerX = cord(1);
centerY = cord(2);
centerZ = cord(3);

radius = r;
circle= (YInImage - centerY).^2 ...
    + (XInImage - centerX).^2 + (ZInImage - centerZ).^2 <= radius.^2;

end







































