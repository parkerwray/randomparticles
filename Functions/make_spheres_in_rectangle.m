function make_spheres_in_rectangle(cord, r, lower_bound, upper_bound)

%{
    This function plots spheres within a rectangle. The rectangle is also
    plotted to show the boundaries of the simulation
%}

figure, 
transparency = 0.2;
plotcube(upper_bound-lower_bound, lower_bound, transparency);
hold on 

[x,y,z] = sphere;
% figure,
hold on 
for i = 1:size(cord,1)
    surf(r(i).*x+cord(i,1), r(i).*y+cord(i,2), r(i).*z+cord(i,3))
    %shading interp
end
title(['Number of particles: ', num2str(size(cord,1))])

end