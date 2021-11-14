function [E] = find_E_fast(X, D)
%This function returns 1 row sparse matrix E that matches for each pixel in X a
%signature from D using matched filter
X = reshape(X, size(X, 1)*size(X, 2), size(X, 3));
X = X.';
%normalize
norms_D = sqrt(sum(abs(D).^2, 2));
norms_D_mat = repmat(norms_D, 1, size(D, 2));
D_normed = D ./ norms_D_mat;
clear norms_D_mat

norms_X = sqrt(sum(abs(X).^2, 1));
norms_X_mat = repmat(norms_X, size(X, 1), 1);
X_normed = (X) ./ norms_X_mat;

clear norms_X_mat
X_normed(isnan(X_normed)) = 0;

ind = 1:length(norms_X);
clear norms_X
E = sparse(size(X, 2), size(D, 1));
val_mats = zeros(1, size(E, 1));

[val_mats(ind), max_ind] = max(real(D_normed*X_normed(:, ind)));
Ind = sub2ind(size(E), ind, max_ind);
E(Ind) = max((real(diag(D_normed(max_ind, :)*X(:, ind))) ./ (norms_D(max_ind))).*(abs(val_mats) > 0.9)', 0); %98
E = real(E);

end
