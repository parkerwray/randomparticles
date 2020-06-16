
function test_PBCs


r = 20;
[x,y,z] = sphere(r);
x = r.*x;
y = r.*y;
z = r.*z;

lower_bounds = [-2,-2,-2].*r;
upper_bounds = -lower_bounds;

sim_size = upper_bounds-lower_bounds;

cord = [sim_size(1)/2, sim_size(2)/2, sim_size(3)/2];
for i = 1:200
    
    cord = cord+[1,1,1];
    cord = periodic_BC_3D(cord, lower_bounds, upper_bounds);
    new_cord = make_mirror_3D(cord, r, lower_bounds, upper_bounds);
    
    clf;
    for j = 1:size(new_cord,1)  
        hold on 

        surf(x+new_cord(j,1), y+new_cord(j,2), z+new_cord(j,3))
        shading interp
        xlim([lower_bounds(1),upper_bounds(1)]);
        ylim([lower_bounds(2),upper_bounds(2)]);
        zlim([lower_bounds(2),upper_bounds(3)]);
        view([-10,-40,30])

    end

    F(i) = getframe;

end
movie(F,500,10)



end

