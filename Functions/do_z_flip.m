function cords = do_z_flip(cords)


top_cords = cords;
bottom_cords = [1,1,-1].*cords;

cords = [top_cords; bottom_cords];
idx_origin_particle = find(cords(:,1)==0& cords(:,2) == 0 & cords(:,3) == 0);
cords(idx_origin_particle,:) = [];
cords = unique(cords,'rows');   
cords = [[0,0,0];cords];


end