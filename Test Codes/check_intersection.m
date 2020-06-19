function flag = check_intersection(cords, radii)
    flag = 0;
    for i = 1:size(cords, 1)
       for j = 1:size(cords, 1)
          if (i ~= j)
             if norm(cords(i,:)-cords(j,:)) < radii(i)+radii(j) 
                dist = norm(cords(i,:)-cords(j,:));
                sum_rad = radii(i) + radii(j);
                flag = 1;
                return;
             end
          end
       end
    end
end