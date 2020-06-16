function [cords, lower_bound, upper_bound] = do_flips(cords,  lower_bound, upper_bound)

v1 = [1,0,0]; %Y-Z plane orth unit vector
v2 = [0,1,0]; %X-Z plane orth unit vector
v3 = [0,0,1]; %X-Y plane orth unit vector

for idx = 1:size(cords, 1)
    cords = [cords; cords(idx,:)-2.*dot(v1,cords(idx,:)).*v1];
    cords = [cords; cords(idx,:)-2.*dot(v2,cords(idx,:)).*v2];
    cords = [cords; cords(idx,:)-2.*dot(v1,cords(idx,:)).*v1-2.*dot(v2,cords(idx,:)).*v2];
    cords = [cords; cords(idx,:)-2.*dot(v3,cords(idx,:)).*v3];
    cords = [cords; cords(idx,:)-2.*dot(v3,cords(idx,:)).*v3-2.*dot(v1,cords(idx,:)).*v1];
    cords = [cords; cords(idx,:)-2.*dot(v3,cords(idx,:)).*v3-2.*dot(v2,cords(idx,:)).*v2];   
    cords = [cords; cords(idx,:)-2.*dot(v3,cords(idx,:)).*v3-2.*dot(v1,cords(idx,:)).*v1-2.*dot(v2,cords(idx,:)).*v2];
end
cords = unique(cords,'rows');    

xmin = min(cords(:,1));
ymin = min(cords(:,2));
zmin = min(cords(:,3));

xmax = max(cords(:,1));
ymax = max(cords(:,2));
zmax = max(cords(:,3));

lower_bound = [xmin, ymin, zmin];
upper_bound = [xmax, ymax, zmax];

idx_remove = [];
for idx = 1:size(cords,1)
    if cords(idx,1) == xmax || cords(idx,2) == ymax || cords(idx,3) == zmax
        idx_remove = [idx_remove, idx];
    end
end
    
cords(idx_remove,:) = [];

end
