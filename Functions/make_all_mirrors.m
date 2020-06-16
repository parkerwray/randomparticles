function final_cords = make_all_mirrors(cords, r, lower_bound, upper_bound)
    
    final_cords = [];    
    for idx2 = 1:size(cords,1)
        dummy = make_mirror_3D(cords(idx2,:), r, lower_bound, upper_bound);
        final_cords = [final_cords;dummy];
    end


end