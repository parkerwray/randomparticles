function [r, final_cords] = make_all_mirrors(cords, radii, lower_bound, upper_bound)
    
    final_cords = [];  
    r = [];
    for idx2 = 1:size(cords,1)
        dummy = make_mirror_3D(cords(idx2,:), radii(idx2), lower_bound, upper_bound);
        final_cords = [final_cords;dummy];
        r = [r;radii(idx2).*ones(size(dummy,1),1)]; %This line was originaly incorrect. -PW
    end


end