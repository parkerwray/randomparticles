%% RUN TESTS
clear;
clc;
% test_PBCs; % Pass
% test_fcc; % Pass
% test_flips; % Pass
% test_make_random_v2; % Pass
dimension = 3;

r = 100;
sigma = 0;
distr = @(~) random('normal', r, sigma);
ff = 0.8;

margin = 0.01;

bounds = [4,4,2]; 
giggles = 100;
tic;


%[radii, cords, bounds, a, am] = ...
%    make_random_fcc_v2(r, radii, ff, bounds, Nspheres, giggles, dimension);
[cords, bounds, a] = make_fcc_3D(r, bounds);
if dimension == 2
    area = (2*a)^2 * bounds(2,1) * bounds(2,2);
else
    area = ((2*a)^3) * bounds(2,1) * bounds(2,2) * bounds(2,3);  
end
disp(num2str(area))
radii = get_radii(area, ff, distr, margin, dimension);
if dimension == 2
    disp(['Created area fill fraction: ', ...
        num2str(100*get_total_volume(radii(1:size(cords,1)),2) / area)]);
else
    disp(['Created fill fraction: ', ...
        num2str(100*get_total_volume(radii(1:size(cords,1)),3) / area)]);
end
disp(['Using ' , num2str(length(radii)), ' particles.']);
%x = get_total_volume(radii, 3)
Nspheres = length(radii);
[radii, cords] = full_randomize(cords, radii, bounds.*a, ...
    giggles, dimension);
%disp("Actual ff (assuming non-uniformity)" + ...
%    num2str(get_total_volume(unique(radii), dimension)/area))
toc;
figure, 
plot_radii(radii); %pass
has_intersections = check_intersection(cords, radii) %pass
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
disp(['X distance: ', num2str((bounds(2,1)-bounds(1,1)).*a)])
disp(['Y distance: ', num2str((bounds(2,2)-bounds(1,2)).*a)])
disp(['Z distance: ', num2str((bounds(2,3)-bounds(1,3)).*a)])