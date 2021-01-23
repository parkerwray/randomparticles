%%

ropt = [749.2249, 565.8669, 736.8252, 712.1183, 530.1508, 502.0276, 530.1508, 791.2003, 773.4692, 478.2372, ...
  669.5459, 530.1508, 798.5562, 502.0276, 783.2299, 699.5273, 478.2371, 530.1508, 468.0189, 893.6108];

wopt = [0.0599, 0.0379, 0.0578, 0.0395, 0.0410,0.0573,0.0409,0.1054,0.0268, 0.0403, 0.0437, 0.0409, ...
    0.1428, 0.0573, 0.0630, 0.0414, 0.0404, 0.0410, 0.0228, 0.0000];






wopt_cutoff = 10^-4;
min_particle_diff = 1; %in nm
[ropt, wopt] = aggregate_similar_particles(ropt, wopt, wopt_cutoff, min_particle_diff);

type = 'film';
dimension = 2;
bounds = [7, 7, 1];
%giggles = 250;
% Sweep parameters (user plural name convention):
center_radiis = ropt; % linspace(r_mean-r_sigma, r_mean+r_sigma, 3); % spacing set to have r = 140nm, to compare with  previous work
% CHOSE CENTER_RADII SPACING CAREFULLY TO MAXIMIZE USEFULLNESS IN STUDYING THE RESULT! 
fill_fraction = 0.443537;

scale_radius = min(ropt);

margin = 0.01;
reps = 25;
%%
clearvars h p
for i = 1:100
[h(i), p(i), Nspheres, ffs] = test_radii_generation(bounds, scale_radius, ...
    ropt, wopt, fill_fraction, margin, dimension, reps);
end
disp(mean(h));