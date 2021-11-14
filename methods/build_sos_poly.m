function mu = build_sos_poly(U, s, q, overres, filter_siz, filter_siz2, ind_filter, ind_filter2)
%Function to build sum-of-squares annihilation weights
normfac = prod(overres) / prod(filter_siz);
mu_small = gpuArray.zeros(filter_siz2);
for j = 1:length(s)
    filter_pad = gpuArray.zeros(overres);
    filter_pad(ind_filter) = ifftshift(reshape(U(:, j), filter_siz));
    filter = reshape(filter_pad(ind_filter2), filter_siz2);
    mu_small = mu_small + ((1 / s(j))^q) * (abs(ifftn(filter)).^2);
end
muhat_small = fftn(mu_small);
muhat = gpuArray.zeros(overres);
muhat(ind_filter2) = muhat_small;
mu = ifftn(muhat) / normfac;
end