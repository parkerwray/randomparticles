function [new_radii, final_cords] = make_all_mirrors(cords, radii, lower_bound, upper_bound)
    
    final_cords = [];    
    new_radii = [];
    for idx2 = 1:size(cords,1)
        dummy = make_mirror_3D(cords(idx2,:), radii(idx2), lower_bound, upper_bound);
        final_cords = [final_cords;dummy];
        for i = 1:size(dummy, 1)
           new_radii = [new_radii; radii(idx2)];
        end
    end


end