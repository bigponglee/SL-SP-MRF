function result = AhAsc_ADMM(x, overres, Samp, ind_full)
%AtA
res_ov = gpuArray.zeros(overres);
na = overres(1);
nb = overres(2);
nc = overres(3);
x = ifft(reshape(x, na, nb, nc), [], 3);
x = x(ind_full); %还原到res维度
result = Samp .* reshape(x, size(Samp));
res_ov(ind_full) = reshape(result, size(x));
res_ov = fft(reshape(res_ov, overres), [], 3);
result = res_ov(:);
end
