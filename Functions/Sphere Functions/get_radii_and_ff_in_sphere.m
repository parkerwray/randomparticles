function [radii, ff_created, Nspheres] = get_radii_and_ff_in_sphere(scale, ...
    r, center_r, ff, distr, margin, dimension, loud)
    %{
        Generates random radii according to the input distribution
        such that the fill fraction is within the margin
    %}

    if dimension == 2
        area = pi * ((scale * r) ^ 2);
    else
        area = (4/3) * pi * ((scale * r) ^ 3);  
    end
    %disp(num2str(area));

    radii = [center_r];
    while get_total_volume(radii, dimension) < ff * area
        radii = [radii distr(0)];
    end
    total_area = get_total_volume(radii, dimension);
    while abs((total_area / area) - ff) > margin
       %disp(num2str(length(radii)))
       if total_area / area > ff
          radii(1+randi(length(radii)-1)) = [];
       else
          radii = [radii distr(0)];
       end
       total_area = get_total_volume(radii, dimension);
    end
    if loud
        disp(['Requested fill fraction: ', num2str(100*ff)]);
    end
    ff_created = total_area/area;
    if loud
        if dimension == 2
            disp(['Created AREA fill fraction: ', ...
                num2str(100*ff_created)]);
        else
            disp(['Created VOLUME fill fraction: ', ...
                num2str(100*ff_created)]);
        end
        disp(['Using ' , num2str(length(radii)), ' particles.']);
    end
    %x = get_total_volume(radii, 3)
    Nspheres = length(radii);
    
    frac_diff = (ff_created-ff)/ff;
    if loud && abs(frac_diff) > margin
        for idx = 1:5
        disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
        disp(['ERROR! The fill fraction you requested was not satisfied!'])
        end
    end
    if loud
        disp(['The fractional difference in', newline,...
            'created and requestion fill fraction is: ',...
            num2str(100.*frac_diff), '%.'])
        disp(['The error margin you requested was: ', num2str(100.*margin),'%'])
    
    
    disp(newline)
    end
end