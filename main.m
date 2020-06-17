%% RUN TESTS
clear;
clc;
% test_PBCs; % Pass
% test_fcc; % Pass
% test_flips; % Pass
% test_make_random_v2; % Pass
dimension = 3;

r = 500;
sigma = 100;
distr = @(~) random('normal', r, sigma);
ff = 0.40;

margin = 0.01;

bounds = [3,3,2]; 
giggles = 100;
tic;

area = (2000^3) * bounds(1) * bounds(2) * bounds(3);
disp(num2str(area))
disp("boohoo")
radii = get_radii(area, ff, distr, margin, 3);
Nspheres = length(radii);
disp("aha")
[radii, cords, bounds, a, am] = ...
    make_random_fcc_v2(r, radii, ff, bounds, Nspheres, giggles, dimension);
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
make_spheres(cords, radii, bounds(1,:), bounds(2,:));
disp(['X distance: ', num2str(bounds(2,1)-bounds(1,1))])
disp(['Y distance: ', num2str(bounds(2,2)-bounds(1,2))])
disp(['Z distance: ', num2str(bounds(2,3)-bounds(1,3))])