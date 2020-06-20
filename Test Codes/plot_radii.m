function [mu, sigma] = plot_radii(radii)

%histogram(unique(radii)); % WHY UNIQUE?????
histogram(radii)
mu = mean(radii);
sigma = std(radii);
title(['Mean: ', num2str(mu), ' Sigma: ', num2str(sigma)])

end