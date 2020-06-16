function [cords, lower_bounds, upper_bounds] = center_cords(cords,r, lower_bounds, upper_bounds)

dummy = make_all_mirrors(cords, r, lower_bounds, upper_bounds);
lower_bounds = lower_bounds-mean(dummy);
upper_bounds = upper_bounds-mean(dummy);
cords = cords-mean(dummy);

[~, test_idx] = min(sum(abs(cords),2));
cords(test_idx,:) = round(cords(test_idx,:));

dummyx = cords(:,1);
dummyx(abs(dummyx)<10^(-12)) = 0;

dummyy = cords(:,2);
dummyy(abs(dummyy)<10^(-12)) = 0;

dummyz = cords(:,3);
dummyz(abs(dummyz)<10^(-12)) = 0;

cords = [dummyx, dummyy, dummyz];
end
