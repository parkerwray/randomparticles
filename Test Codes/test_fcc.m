function [ff_calc, Ntouch, ff, is_zero] = test_fcc

clear;
clc;
ff=0.1;

ff = 0.1:0.01:0.73;

parfor (idx = 1:length(ff))
    lower_bound = [0,0,0];
    upper_bound = [4,4,4];
    r = 10;
    
    [cords, bounds, ~, ~] = make_fcc_v2(r, ff(idx), lower_bound, upper_bound);
    [cords, lower_bound, upper_bound] = center_cords(cords, r, bounds(1,:), bounds(2,:));
    cord_final = make_all_mirrors(cords, r, lower_bound, upper_bound);

    dummy = cord_final;
    norm_dist = check_distance_function(dummy, r, lower_bound, upper_bound);
    values = norm_dist(norm_dist<1);
    values(values == 0) = [];
    Ntouch(idx) = sum(values<1);
    [~, ff_calc(idx)] = ...
        visualize_spheres(dummy,...
        round(upper_bound(1)), round(upper_bound(3)),...
        r, ff(idx), 0);
    dist = sqrt(dummy(:,1).^2+dummy(:,2).^2+dummy(:,3).^2);
    
    is_zero(idx,:) = min(abs(dist));
    
    
    
end

figure, 
plot(ff, ff_calc)
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

end