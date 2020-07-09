%% RUN TESTS
clear;
clc;
% test_PBCs; % Pass
% test_fcc; % Pass
% test_flips; % Pass
% test_make_random_v2; % Pass
loud = 1;
dimension = 2;
type = "film";
scale = 10;
r = 100;
sigma = 10;
distr = @(~) random('normal', r, sigma);
ff = 0.4;

center_radius = 100;

margin = 0.01;

bounds = [2,2,2]; 
giggles = 100;
tic;
if strcmp(type, "sphere") == 1
    [radii, ff, Nspheres] = get_radii_and_ff_in_sphere(scale, r, ...
        center_radius, ff, distr, margin, dimension, loud);
    cords = full_randomize_in_sphere(radii, scale*r, giggles, dimension, loud);
elseif strcmp(type, "film") == 1
    if dimension == 3
        [cords, bounds, a] = make_fcc_3D(r, bounds, loud);
    else
        [cords, bounds, a] = make_fcc_2D(r, bounds, loud);
    end
    [radii, ff, Nspheres] = get_radii_and_ff(bounds, a,...
        center_radius, ff, distr, margin, dimension, loud);
    if dimension == 2
       max_radii = max(radii(:));
       bounds(:,3) = bounds(:,3).*max_radii/a/2;
    end
    [radii, cords] = full_randomize(cords, radii, bounds.*a, ...
        giggles, dimension, loud);
else
    disp("Invalid geometry requested."); 
end
if loud
    toc;
    figure, 
    plot_radii(radii); %pass
    has_intersections = check_intersection(cords, radii); %pass
end
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
if loud
    
if strcmp(type, "film") == 1
    make_spheres(cords, radii, bounds(1,:).*a, bounds(2,:).*a);
    disp('Simulation region:')
    disp(['X distance: ', num2str((bounds(2,1)-bounds(1,1)).*a)])
    disp(['Y distance: ', num2str((bounds(2,2)-bounds(1,2)).*a)])
    disp(['Z distance: ', num2str((bounds(2,3)-bounds(1,3)).*a)])
elseif strcmp(type, "sphere") == 1
    make_spheres_in_sphere(cords, radii, r*scale);
    %make_spheres(cords, radii, [-1000,-1000,-1000], [1000,1000,1000]);
end

end