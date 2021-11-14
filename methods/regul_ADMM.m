function res = regul_ADMM(x, MtM, overres, gam, lambda)
x = reshape(x, overres);
temp1 = MtM .* x;
res = gam * lambda * temp1(:);
clear temp1
end
