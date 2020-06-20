%% RUN TESTS
clear;
clc;
% test_PBCs; % Pass
% test_fcc; % Pass
% test_flips; % Pass
% test_make_random_v2; % Pass
dimension = 3;

r = 100;
sigma = 10;
distr = @(~) random('normal', r, sigma);
ff = 0.4;

margin = 0.01;

bounds = [2,2,2]; 
giggles = 100;
tic;

[cords, bounds, a] = make_fcc_3D(r, bounds);
[radii, ff, Nspheres] = get_radii_and_ff(bounds, a,...
    ff, distr, margin, dimension);
[radii, cords] = full_randomize(cords, radii, bounds.*a, ...
    giggles, dimension);

toc;
figure, 
plot_radii(radii); %pass
has_intersections = check_intersection(cords, radii); %pass
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
make_spheres(cords, radii, bounds(1,:).*a, bounds(2,:).*a);
disp('Simulation region:')
disp(['X distance: ', num2str((bounds(2,1)-bounds(1,1)).*a)])
disp(['Y distance: ', num2str((bounds(2,2)-bounds(1,2)).*a)])
disp(['Z distance: ', num2str((bounds(2,3)-bounds(1,3)).*a)])