

dimension = 2;
boundary_type = "sphere";
r = [50,100,150,200];
probability = [0.25, 0.25, 0.25, 0.25];
ff = 0.4;
center_radius = 100;
size = 10000;
giggles = 100;
loud = 1;
sim_folder = "/home/parkerwray/hypnos/Codes/randomparticles";
save_folder = "/home/parkerwray/hypnos/Codes/randomparticles";

[spheres, ff] = make_random_particle_film(dimension, boundary_type, r, probability, ff, center_radius, size, giggles, loud, sim_folder, save_folder);










