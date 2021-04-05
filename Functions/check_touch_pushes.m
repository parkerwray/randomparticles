function pushes = check_touch_pushes(new_cord, cords, r, radii, lower_bounds, upper_bounds)

% This function determines if a particle (centered at new_cord, with radius
% r) is touching any other particles (given by particle list cords, and
% radius r)

pushes = zeros(size(cords, 1), 3);
for idx = 1:size(cords,1)
    %dist = get_distance(new_cord, cords(idx, :), lower_bounds, upper_bounds);
    %dist = norm(new_cord - cords(idx,:));
    [dist, dcord] = get_distance(new_cord, cords(idx,:), lower_bounds, upper_bounds);
    if dist < r + radii(idx) && dist > 0% - (10^(-3))
        %disp(dist);
        %pushes(idx,:) = new_cord - cords(idx, :);
        pushes(idx,:) = dcord;
        %pushes(idx,:) = pushes(idx,:) ./ norm(pushes(idx,:));
        %pushes(idx,:) = pushes(idx,:)/10 + ((pushes(idx,:)/norm(pushes(idx,:))));
        %pushes(idx,:) = (pushes(idx,:).*(r+radii(idx)))./(norm(pushes(idx,:))^2);
        %    ((pushes(idx,:)./norm(pushes(idx,:))));
        pushes(idx,:) = pushes(idx,:)./norm(pushes(idx,:)).*(r+radii(idx)-dist+0.1);
    end
end

end