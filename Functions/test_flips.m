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