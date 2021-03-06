


clc;
bounds = [2,2,2]; %MUST BE AN INTERGER > 1. 
% Note: bounds = [1,1,1]; Should work. But I find not all particles are
% generated. This may be because the for loop?? Also, at this point, bounds
% only gives 1/4th of the simulation region. The code flips this particle
% diftribution across x, y, and z such that the center particle is always
% in the middle of the entire simulation region (this was a solution to an
% original headache...)

% Interger and interger + 1/2 bounds will create a proper periodic
% structure. BUT, for the EM calculations, it is necessary to have a
% particle always at the origin. This is our "observation particle". 
% If you use interger bounds this origin particle is always properly
% alligned with the lattice. If you use interger + 0.5 you are not
% alligned. This can be easily seen by using dimension = 2.

% NOTE: The necessity that the origin particle never move is important. 
% I believe we talked about this when we first met, but I forgot to
% re-mention it again when we picked back up the project. I apologize for
% that. I make the code such that the origin particle is always the first 
% row element. This is easy to skip it when moving particles. You should
% see this in the random move funciton. 

dimension = 3;
r = 10;

[cord, bounds, a, Nspheres] = ...
    make_fcc_v3(r, bounds, dimension);

r2 = r.*ones(size(cord,1),1); %account for all same radius using updated functions

% % This check distance function does not have the capability to deal with 
% % different radii. I'm unsure if you have updated this already.  
% [overlap_distance, overlap_idx] = ...
%     check_distance_function_v2(cord, r, bounds(1,:).*a, bounds(2,:).*a);


[r3, final_cord] =...
    make_all_mirrors(cord, r2, bounds(1,:).*a, bounds(2,:).*a);

make_spheres_in_rectangle(final_cord, r3, bounds(1,:).*a, bounds(2,:).*a);

[overlap_distance, overlap_idx] = ...
    check_distance_function_v2(final_cord, r, bounds(1,:).*a, bounds(2,:).*a);

% NOTE: Nspheres is not based on final_cord. This is because final_cord has
% the coordinates of the periodic repeated particles at the x-y edges. These 
% coordinates are necessary for plotting and EM simulations. But the total 
% number of physical particles is not changing. 
Vsphere = (4/3)*pi*(r^3);
Vspace = prod((bounds(2,:)-bounds(1,:)).*a);
ff = (Nspheres).*Vsphere./Vspace;

% FF is slightly off from theoretical limit based on NO Z peridicity. 
