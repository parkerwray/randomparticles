function delta = get_all_move_directions_sphere(cords, radii, R)
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
    
           pushes = check_touch_pushes_sphere(cords(idx1,:), dummy, ...
               radii(idx1), dummyr, R);
           for idx = 1:size(pushes, 1)
              delta(idx1,:) = delta(idx1,:) + pushes(idx,:); 
           end
       end
    end
end