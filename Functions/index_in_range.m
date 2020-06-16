function mat2 = index_in_range(mat, lda_set, loud)

lda = mat.lda;
if min(lda_set)< min(mat.lda)
    disp('Error: Lowerbound wavelength is not supported by data of this material!')
end
if max(lda_set)> max(mat.lda)
    disp('Error: Upperbound wavelength is not supported by data of this material!')
end

mat2.n = interp1(mat.lda, mat.n, lda_set);
mat2.k = interp1(mat.lda, mat.k, lda_set);
mat2.e1 = interp1(mat.lda, mat.e1, lda_set);
mat2.e2 = interp1(mat.lda, mat.e2, lda_set);


mat2.e = mat2.e1+1i.*mat2.e2;
mat2.m = mat2.n+1i.*mat2.k;
mat2.lda = lda_set;

if loud
    [~, idx_min] = min(abs(mat.lda-min(lda_set)));
    [~, idx_max] = min(abs(mat.lda-max(lda_set)));
    figure, 
    subplot(1,2,1)
    hold on 
    plot(mat.lda(idx_min:idx_max), mat.n(idx_min:idx_max))
    plot(mat.lda(idx_min:idx_max), mat.k(idx_min:idx_max))
    hold off
    xlabel('Wavelength')
    ylabel('index')
    legend('n','k');
    subplot(1,2,2)
    hold on 
    plot(mat2.lda, mat2.n)
    plot(mat2.lda, mat2.k)
    hold off
    xlabel('Wavelength')
    ylabel('Interpolated index')
    legend('n','k');
end





end










