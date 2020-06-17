function total = get_total_volume(radii, dimension)
    %{
        Gets the total volume/area covered by the given
        set of radii
    %}
    total = 0;
    for i = 1:length(radii)
        if (dimension == 2)
           total = total + pi * (radii(i) ^ 2);
        else
            total = total + pi * (4 / 3) * (radii(i) ^ 3);
        end
    end
end