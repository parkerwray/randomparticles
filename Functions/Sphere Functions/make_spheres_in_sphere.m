function make_spheres_in_sphere(cord, r, R)

% Generates images for spheres

figure, 
transparency = 0.2;
%{
[x,y,z] = sphere(R);
h = surfl(x, y, z); 
set(h, 'FaceAlpha', transparency)
hold on
%}
plotcube([2*R,2*R,2*R], [-R,-R,-R], 0);
hold on

[x,y,z] = sphere;
% figure,
% hold on 
for i = 1:size(cord,1)
    surf(r(i).*x+cord(i,1), r(i).*y+cord(i,2), r(i).*z+cord(i,3))
    %shading interp
end
title(['Number of particles: ', num2str(size(cord,1))])
hold off
end