
%% Make Randomly positioned same sized spheres inside of a rectangle
r = 40;
ff = 0.3;
bounds = [100, 100, 50]./r;
dimension = 3;
giggles = 1;

[cord, bounds, a, am, Nspheres] = ...
    make_random_fcc_v2(r, ff, bounds, giggles, dimension);
lower_bound = bounds(1,:);
upper_bound = bounds(2,:);
make_spheres_in_rectangle(cord, r, lower_bound, upper_bound)

%% Get material refractive index for spheres
mat_file = 'Z:\Project - Dusty Plasma MURI\Material Data\Refractive Index Info\Ag_Yang.csv';
name = 'mat';
mat = pw_read_refinfo_mat(name, mat_file, 1);

%% Compile sphere [radius, x, y, z, n, k] into matrix 
wavelengths = 300:1:400;
dummy = [r.*ones(size(cord,1),1), cord];
mat = index_in_range(mat, wavelengths, 1);

for idx = 1:10 
  sphere_data(:,:,idx) = [dummy,...
        mat.n(idx).*ones(size(cord,1),1),...
        mat.k(idx).*ones(size(cord,1),1)];
end