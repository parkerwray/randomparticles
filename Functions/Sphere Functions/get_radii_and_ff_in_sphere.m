function [radii, ff_created, Nspheres] = get_radii_and_ff_in_sphere(scale, ...
    r, ff, distr, margin, dimension)
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

    radii = [];
    while get_total_volume(radii, dimension) < ff * area
        radii = [radii distr(0)];
    end
    total_area = get_total_volume(radii, dimension);
    while abs((total_area / area) - ff) > margin
       %disp(num2str(length(radii)))
       if total_area / area > ff
          radii(randi(length(radii))) = [];
       else
          radii = [radii distr(0)];
       end
       total_area = get_total_volume(radii, dimension);
    end
    disp(['Requested fill fraction: ', num2str(100*ff)]);
    
    ff_created = total_area/area;
    if dimension == 2
        disp(['Created AREA fill fraction: ', ...
            num2str(100*ff_created)]);
    else
        disp(['Created VOLUME fill fraction: ', ...
            num2str(100*ff_created)]);
    end
    disp(['Using ' , num2str(length(radii)), ' particles.']);
    %x = get_total_volume(radii, 3)
    Nspheres = length(radii);
    
    frac_diff = (ff_created-ff)/ff;
    if abs(frac_diff) > margin
        disp(['ERROR! The fill fraction you requested was not satisfied!'])
    end
    disp(['The fractional difference in', newline,...
        'created and requestion fill fraction is: ',...
        num2str(100.*frac_diff), '%.'])
    disp(['The error margin you requested was: ', num2str(100.*margin),'%'])
    
    
    
    disp(newline)
end