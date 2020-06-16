function [dist, dcord] = get_distance(cord1, cord2, lower_bounds, upper_bounds)

dcord1 = periodic_BC_3D((cord1-cord2), lower_bounds, upper_bounds);
dist1 = sqrt(dcord1(1)^2+dcord1(2)^2+dcord1(3)^2);

dcord2 = periodic_BC_3D((cord2-cord1), lower_bounds, upper_bounds);
dist2 = sqrt(dcord2(1)^2+dcord2(2)^2+dcord2(3)^2);

if dist2 < dist1
    dist = dist2;
    dcord = dcord2;
else
    dist = dist1;
    dcord = dcord1;
end


end
