function cords = fix_overlap(cords, radii, lower_bounds, upper_bounds)
    %{ 
        Eliminates overlap caused by changing radii
    %}

    disp(['Moving particles away from each other to account for overlap',...
    newline, 'caused by changing the particle radii.'])
    disp('Iteration number:')
    lineLength = fprintf(num2str(1),21.1);


    count = 0;
    delta = get_all_move_directions(cords, radii, ...
        lower_bounds, upper_bounds);
    while ~isequal(delta, zeros(length(cords), 3))
        count = count + 1;
        if mod(count, 10) == 0
            fprintf(repmat('\b',1,lineLength))
            lineLength = fprintf(num2str(count),21.1);
            pause(0.001);
%             disp(count)
        end
        %disp("moving!")
        cords = cords + delta;
        % Fix periodicity/outofbounds
        for idx = 1:size(cords, 1)
           
            new_cord0 = periodic_BC_3D(cords(idx,:), lower_bounds + ...
                [0,0,radii(idx)], upper_bounds - [0,0,radii(idx)]);
            %new_cord = make_mirror_3D(new_cord0, r(idx), ...
            %    lower_bounds, upper_bounds); 
            
            cords(idx, :) = new_cord0;
        end
        delta = get_all_move_directions(cords, radii, ...
            lower_bounds, upper_bounds);
    end
    disp(newline)
end