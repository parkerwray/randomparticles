function radii = get_radii(area, ff, distr, margin, dimension)
    %{
        Generates random radii according to the input distribution
        such that the fill fraction is within the margin
    %}

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
    if dimension == 2
        disp(['Created area fill fraciton: ', ...
            num2str(100*total_area / area)]);
    else
        disp(['Created fill fraction: ', ...
            num2str(100*total_area / area)]);
    end
    disp(['Using ' , num2str(length(radii)), ' particles.']);
end