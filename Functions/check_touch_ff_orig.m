


function [Ntouch, ffcalc, is_zero, dist, overlap_idx] = check_touch_ff_orig(cord, ff, r, lower_bound, upper_bound)

    dummy = cord;
    [overlap_distance, overlap_idx] = check_distance_function_v2(dummy, r, lower_bound, upper_bound);
    Ntouch = size(overlap_idx,1);
%     [I, ffcalc] = ...
%         visualize_spheres(dummy,...
%         round(upper_bound(1)-lower_bound(1)), round(upper_bound(3)-lower_bound(3)),...
%         r, ff, 0);
%     figure, 
%     imshow(squeeze(I(:,:,5)))
    
    ffcalc = ff;
    dist = sqrt(dummy(:,1).^2+dummy(:,2).^2+dummy(:,3).^2);
    
    is_zero = min(abs(dist));

    disp(['Desired fill fraction: ', num2str(100*ff)])
    disp(['Calculated fill fraction: ', num2str(100*ffcalc)])
    disp(['Number of touching spheres: ', num2str(Ntouch)])
    
    if Ntouch > 0
        disp(['Worst fractional overlap: ', num2str(100.*min(overlap_distance)), '%'])
        disp(['Overlaped spheres: '])
            disp(num2str(overlap_idx))
            keyboard;
    end
    
    disp(['Orig particle distance to orig: ', num2str(is_zero)])
    disp(newline)
end