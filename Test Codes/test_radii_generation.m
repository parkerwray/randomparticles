function [h, p, Nspheres, ffs] = test_radii_generation(bounds, scale_radius, ...
    radii, weights, ff, margin, dimension, reps)

%typical margin of 0.01
a = (4 / sqrt(2)) * scale_radius;
bounds = [-bounds; bounds];


Nspheres = zeros(reps * length(radii), 1);
ffs = zeros(reps * length(radii), 1);




single_particle_weights = weights ./ (radii .^ 2);
single_particle_weights = single_particle_weights / sum(single_particle_weights);
dist = @(~) randsample(radii, 1, true, single_particle_weights);

cdf = zeros(length(radii), 2);
cdf(:,1) = radii;
for i = 1:length(single_particle_weights)
    cdf(i,2) = sum(single_particle_weights(1:i));
end


all_radii = [];


for i = 1:length(radii)
    for j = 1:reps
        [these_radii, ffs((i - 1) * reps + j), Nspheres((i-1) * reps + j)] = ...
            get_radii_and_ff(bounds, a, radii(i), ff, dist, margin, dimension, 0);
        all_radii = [all_radii, these_radii(2:end)];
    end
end
%{
figure,
histogram(Nspheres);
figure,
histogram(ffs);
%}
%{
figure,
hold on
histogram(all_radii, 10*length(radii), 'Normalization', 'probability');
plot(radii, single_particle_weights);
hold off
%}
[h, p] = chi2gof(all_radii, 'Ctrs', cdf(:,1), 'CDF', {@cdf_eval, cdf});

end

function cdf = cdf_eval(x, cdf_in)
    idx = zeros(length(x), 1);
    for i = 1:length(x)
        idx(i) = find(cdf_in(:,1) - x(i) >= 0, 1);
        if x(i) ~= cdf_in(idx(i), 1)
            idx(i) = idx(i) - 1;
        end
    end
    cdf = cdf_in(idx, 2);
    %disp(cdf)
end