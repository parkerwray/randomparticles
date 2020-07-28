function FLAG = check_fill_fraction(bounds, radii, ff, margin, dimension)
    FLAG = 0;

    if dimension == 2
        area = (bounds(2,1)-bounds(1,1)) * (bounds(2,2)-bounds(1,2));
    else
        area = (bounds(2,1)-bounds(1,1)) * (bounds(2,2)-bounds(1,2)) * (bounds(2,3)-bounds(1,3));  
    end
    
    total_area = get_total_volume(radii, dimension);
    
    if abs((total_area / area) - ff) > margin && ff > 0
        FLAG = 1;
        disp("Incorrect fill fraction");
        keyboard;
    end
end