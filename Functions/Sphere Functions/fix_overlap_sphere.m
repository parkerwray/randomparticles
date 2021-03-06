function cords = fix_overlap_sphere(cords, radii, R, loud)
    %{ 
        Eliminates overlap caused by changing radii
    %}
    if loud
        disp(['Moving particles away from each other to account for overlap',...
        newline, 'caused by changing the particle radii.'])
        disp('Iteration number:')
        lineLength = fprintf(num2str(1),21.1);
    end

    count = 0;
    delta = get_all_move_directions_sphere(cords, radii, R);
    while ~isequal(delta, zeros(length(cords), 3))
        count = count + 1;
        if loud
            if mod(count, 10) == 0
                fprintf(repmat('\b',1,lineLength))
                lineLength = fprintf(num2str(count),21.1);
                pause(0.001);
    %             disp(count)
            end
        end
        %disp("moving!")
        cords = cords + delta;
        delta = get_all_move_directions_sphere(cords, radii, R);
    end
    if loud
    disp(newline)
    end
end