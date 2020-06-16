function [dist2, overlap_idx] = check_distance_function(cords, r, lower_bounds, upper_bounds)
overlap_idx = [];
for idx = 1:size(cords,1)
    for idx2 = 1:size(cords,1)
        dist(idx, idx2) = get_distance(cords(idx2, :), cords(idx, :), lower_bounds, upper_bounds);
        dist2(idx,idx2) = dist(idx,idx2)./(2*r);
        if dist2(idx,idx2) < 1 && dist2(idx, idx2) > 0
            overlap_idx = [overlap_idx;[idx,idx2]];
        end
    end
end



end