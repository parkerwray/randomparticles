function repeat = check_repeat(sphere_cords)

    last = length(sphere_cords);
    repeat = 0;
    for idx=1:last-1
        if isequal(sphere_cords{idx}, sphere_cords{last})
            repeat = 1;
            break;
        end
    end
    
end