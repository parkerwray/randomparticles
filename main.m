%% RUN TESTS
clear;
clc;
% test_PBCs; % Pass
% test_fcc; % Pass
% test_flips; % Pass
% test_make_random_v2; % Pass
dimension = 2;

r = 500;
sigma = 100;
distr = @(~) random('normal', r, sigma);
ff = 0.4;

margin = 0.01;

bounds = [20,20,1]; 
giggles = 100;
tic;


%[radii, cords, bounds, a, am] = ...
%    make_random_fcc_v2(r, radii, ff, bounds, Nspheres, giggles, dimension);
[cords, bounds, a] = make_fcc_v3(r, bounds, dimension);
if dimension == 2
    area = (2*a)^2 * bounds(2,1) * bounds(2,2);
else
    area = ((2*a)^3) * bounds(2,1) * bounds(2,2) * bounds(2,3);  
end
disp(num2str(area))
radii = get_radii(area, ff, distr, margin, dimension);
%x = get_total_volume(radii, 3)
Nspheres = length(radii);
disp("a");
[radii, cords] = full_randomize(cords, radii, bounds.*a, ...
    giggles, dimension);
%disp("Actual ff" + num2str(get_total_volume(radii, 3)/area))
toc;
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