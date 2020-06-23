function pushes = check_touch_pushes_sphere(new_cord, cords, r, radii, R)

% This function determines if a particle (centered at new_cord, with radius
% r) is touching any other particles (given by particle list cords, and
% radius r)

pushes = zeros(size(cords, 1)+1, 3);
for idx = 1:size(cords,1)
    dist = norm(new_cord - cords(idx,:));
    if dist < r + radii(idx)
        pushes(idx,:) = new_cord - cords(idx, :);
        %pushes(idx,:) = pushes(idx,:) ./ norm(pushes(idx,:));
        %pushes(idx,:) = pushes(idx,:)/10 + ((pushes(idx,:)/norm(pushes(idx,:))));
        %pushes(idx,:) = (pushes(idx,:).*(r+radii(idx)))./(norm(pushes(idx,:))^2);
        %    ((pushes(idx,:)./norm(pushes(idx,:))));
        pushes(idx,:) = pushes(idx,:)./norm(pushes(idx,:)).*(r+radii(idx)-dist+0.1);
    end
end
if norm(new_cord) + r > R
    pushes(size(pushes,1),:) = (new_cord / norm(new_cord)) * ...
        (R - (norm(new_cord)+r));
end

end