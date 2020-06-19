function new_cord = random_move(cord, r, dimension)
% Randomly move the particle +1, 0, or -1 in the x, y, and z dimensions.
r = floor(r);
if r == 0
    r = 1;
end
if dimension == 2
       new_cord = [cord(1)+randi(r).*round(2*rand-1),...
       cord(2)+randi(r).*round(2*rand-1),...
       cord(3)];
else
   new_cord = [cord(1)+randi(r).*round(2*rand-1),...
       cord(2)+randi(r).*round(2*rand-1),...
       cord(3)+randi(r).*round(2*rand-1)];
end
  
end