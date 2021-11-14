function G = build_gram_matrix(gradx, filter_siz, filter_siz2, ind_filter2)
% Function to build gram matrix G=T(x)^*T(x)
% using ffts and neighborhood operators
n = size(gradx, 4);
for j = 1:n
    gradx_ifft(:, :, :, j) = ifftn(gradx(:, :, :, j));
end
clear gradx

gradx_ifft = gather(gradx_ifft);
g = fftn(sum(conj(gradx_ifft).*gradx_ifft, 4));

g = fftshift(reshape(g(ind_filter2), filter_siz2));
G = im2colstep(real(g), filter_siz, [1, 1, 1]) + 1i * im2colstep(imag(g), filter_siz, [1, 1, 1]);

G = rot90(G, -1);
end