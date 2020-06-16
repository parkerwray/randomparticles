function [sphere_m, medium_m, k, lda] = get_same_index(mat, set_lda, Nspheres, loud)

lda_min = ceil(min(mat.lda));
lda_max = floor(max(mat.lda));

lda = lda_min:1:lda_max;
n = interp1(mat.lda, mat.n, lda);
k = interp1(mat.lda, mat.k, lda);


idx = find(lda==set_lda);
lda = lda(idx); %EACH WAVELENGTH IS A NEW SIMULATION!
sphere_m = n(idx)+1i.*k(idx);  
sphere_m = repmat(sphere_m,Nspheres,1);
medium_m = 1+1i*0;
k = 2.*pi./lda;

if loud == 1
    figure, 
    plot(mat.lda, mat.n, mat.lda, mat.k)
    hold on 
    plot(lda, n, lda, k)
    xlabel('Wavelength (nm)')
end




end
