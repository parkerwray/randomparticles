function test_make_random_v2

clear;
clc;
close all;

ff = 0.1:0.01:0.7;

parfor idx = 1:length(ff)
    lower_bound = [0,0,0];
    upper_bound = [3,3,3];
    r = 10;

    [cord, bounds, ~, ~] = make_fcc_v2(r, ff(idx), lower_bound, upper_bound);
    [cord_flip, lower_bound, upper_bound] = do_flips(cord, bounds(1,:), bounds(2,:));
    
    cord_final_per = make_all_mirrors(cord_flip, r, lower_bound, upper_bound);
    
    cord_random = make_random(cord_flip, r, lower_bound, upper_bound);
    cord_final = make_all_mirrors(cord_random, r, lower_bound, upper_bound);

    [Ntouch(idx), ffcalc(idx), is_zero(idx)] = ...
        check_touch_ff_orig(cord_final_per, ff(idx), r, lower_bound, upper_bound);

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

% [x,y,z] = sphere(r);


% figure,
% hold on 
% for i = 25:26 %size(cord_random,1)
%     surf(r.*x+cord_random(i,1), r.*y+cord_random(i,2), r.*z+cord_random(i,3))
%     shading interp
% end
% title(['Number of particles: ', num2str(size(cord_random,1))])

% figure,
% hold on 
% for i = 1:size(cord_final,1)
%     surf(r.*x+cord_final(i,1), r.*y+cord_final(i,2), r.*z+cord_final(i,3))
%     shading interp
% end
% title(['Number of particles: ', num2str(size(cord_final,1))])

% for idx = 1:size(overlap_idx,1)
%     disp(['Particle 1: ', num2str(cord_final(overlap_idx(idx,1),:))])
%     disp(['Particle 2: ', num2str(cord_final(overlap_idx(idx,2),:))])
%     disp(num2str(norm_dist(overlap_idx(idx,1),overlap_idx(idx,2))))
% 
% end


end
