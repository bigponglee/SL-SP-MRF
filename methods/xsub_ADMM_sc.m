function x_pad = xsub_ADMM_sc(x0, Ahb_pad, mu, mask, Mt, M, MtM, overres, gam, lambda, ind_full, iter, tol)
% ADMM ITER
x = fft(x0, [], 3);
Mx = M(x);
q = gpuArray.zeros(size(M(x)));
Y = gpuArray.zeros(size(M(x)));
ndz = size(q, 4);
resvec = gpuArray.zeros(1, iter);
for i = 1:iter
    % Y subprob
    Z = gam * (Mx + q);
    muinv = repmat((mu + gam).^(-1), [1, 1, 1, ndz]);
    for j = 1:ndz
        Y(:, :, :, j) = fftn(muinv(:, :, :, j).*ifftn(Z(:, :, :, j))); % Y:fft(y)
    end
    % x subprob
    f = @(x)(regul_ADMM(x, MtM, overres, gam, lambda) + AhAsc_ADMM(x, overres, mask, ind_full)); %%%
    rhs = Ahb_pad + lambda * gam * Mt(Y-q);
    % tic;
    [x, ~, ~, ~, ~] = pcg(f, rhs(:), 1e-8, 50, [], [], x(:));

    % L update
    Mx = M(reshape(x, overres));
    res2 = Mx - Y;
    q = q + res2;

    resvec(i) = norm(res2(:)) / norm(Y(:));
    if (iter > 2) && (resvec(i) < tol)
        x_pad = ifft(reshape(x, overres), [], 3);
        return;
    end
end
x_pad = ifft(reshape(x, overres), [], 3);
end
