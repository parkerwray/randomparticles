function delta = get_all_move_directions(cords, radii, lower_bounds, upper_bounds)
    %{
      Finds move directions for particles to reduce overlap 
    %}

    delta = zeros(length(cords),3);
    for idx1 = 1:size(cords, 1)
       if ~isequal(cords(idx1,:), [0,0,0])
           dummy = cords;
           dummy(idx1,:) = [];
           dummyr = radii;
           dummyr(idx1) = [];


            new_cord = make_mirror_3D(cords(idx1,:), radii(idx1), ...
                lower_bounds, upper_bounds);

           for idx2 = 1:size(new_cord, 1)
              pushes = check_touch_pushes(new_cord(idx2,:), dummy, ...
                  radii(idx1), dummyr, lower_bounds, upper_bounds);
              for idx = 1:size(pushes, 1)
                 delta(idx1,:) = delta(idx1,:) + pushes(idx,:); 
              end
           end
       end
       %if norm(delta(idx1,:)) ~= 0
       %   delta(idx1,:) = delta(idx1,:) ./ norm(delta(idx1,:)); 
       %end
    end
end