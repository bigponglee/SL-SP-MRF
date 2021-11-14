function [Xr, Y, m_abs, runtime] = SL_SP(X0, mask, noise, D, noisy)

%% ����˵����
% ���룺X0��ʱ��
% mask:Ƶ�����ģ��  δshift �Ľǵ�Ƶ���ĸ�Ƶ
% noise:Ƶ������
% D���ֵ�
% �����Xr���ؽ����ʱ��  Y:Ƶ���һ�������� m_abs:��һ������

%% Ƶ�򼰹�һ��
m_abs = max(abs(X0(:))); %��һ������
X0 = X0 / m_abs; %��һ��

Y = gpuArray(fft2(X0)); %Ƶ�� ��֡fft
X0 = gather(ifft2(Y));

%% ������
if noisy == 1
    Y = Y + gpuArray(noise);
    disp('Adding noise to the data in Fourier domain');
end

%% ��������
filter_siz = [15, 15, 7]; %�˲�����С
res = size(X0); %����ߴ� N*N*L
filter_siz2 = 2 * filter_siz - [1, 1, 1]; %squared filter dimensions
overres = res + 2 * filter_siz; %over-resolved reconstruction grid
[ind_full, ind_filter, ind_filter2, k] = get_kspaceindices(overres, res, filter_siz); %intialize index sets

lambda = 1e12; %lambda
eta = 10; %epsilon decrease factor;
epsmin = 1e-3; %minimum possible epsilon value
p = 0; % p for schatten p-norm
q = 1 - (p / 2); %q

ADMM_iter = 10; %admm����������
ADMM_tol = 1e-8; %admm�����С���exit tolerance for ADMM
delta = 100; %default ADMM conditioning parameter

maxiter = 10; %giraf����������

%% Ƿ����
Samp = gpuArray(mask); %����ģ��
ind = find(Samp(:)); %��������
[A, At] = defAAt(ind, res);
b = A(Y); %Ƶ��Ƿ����
Ahb = At(b); %��������0  Ƶ��
x_pad = gpuArray.zeros(overres);
x_pad(ind_full) = Ahb; %x
Ahb_pad = fft(x_pad, [], 3); %����������һ�׳���  ��ʱ��ά�Ͻ���fft�任

%% ��������
alpha = res(3) * 1e-2; %%% Treat this as an optimization parameter
dz(:, :, :, 1) = reshape((1i * 2 * pi * (k(1, :))).', overres) / res(2); %%% why divide?- To make it resolution independent
dz(:, :, :, 2) = reshape((1i * 2 * pi * (k(2, :))).', overres) / res(1);
dz(:, :, :, 3) = reshape((1i * 2 * pi * (k(3, :))).', overres) / (alpha);

M = @(z) repmat(z, [1, 1, 1, size(dz, 4)]) .* dz; %M
Mt = @(Z) sum(Z.*conj(dz), 4); %Mt
MtMmask = Mt(M(ones(overres))); %MtM

%% �ֵ��
base = orth(D.');
base = base.';
pinv_D = pinv(base); %Moore-Penrose Pseudoinverse ��
pinv_DD = pinv_D * base;

%%
% figure('Name', 'iter cost', 'NumberTitle', 'off');
cost = zeros(maxiter,1);
SNR = 1:maxiter;
tic; %�����ʱ
for i = 1:maxiter
    xtemp = fft(x_pad, [], 3); %��ʱ��ά�Ͻ���fft�任
    % step 1: Compute sos annihilating polynomial
    gradx = M(xtemp); %Mx
    G = build_gram_matrix(gradx, filter_siz, filter_siz2, ind_filter2); %build gram matrix G=T(x)^*T(x)
    [U, S] = eig(gpuArray(G)); %SVD
    ev = abs(diag(S));
    if i == 1 %initialze epsilon
        eps = 0.001 * max(ev); %auto-init eps
    end
    mu = build_sos_poly(U, ev+eps, q, overres, filter_siz, filter_siz2, ind_filter, ind_filter2); %build sum-of-squares annihilation weights h
    % step 2: ADMM solution of least squares problem
    gama = max(mu(:)) / delta; %gama
    x_pad = xsub_ADMM_sc(x_pad, Ahb_pad, mu, Samp, Mt, M, MtMmask, overres, gama, lambda, ind_full, ADMM_iter, ADMM_tol);
    % step3: projection
    r = ifft2((reshape(x_pad(ind_full), res)));
    r = reshape(r, res(1)*res(2), res(3));
    r = r * pinv_DD;
    x_pad(ind_full) = reshape(fft2(reshape(r, res)), 1, res(1)*res(2)*res(3));
    x = reshape(x_pad(ind_full), res); %ȥ�����
    Xr = gather(ifft2(x)); %ʱ���ؽ�����
    % update epsilon
    eps = max(eps/eta, epsmin);
    % check stopping condition
    cost(i,1) = get_cost(p, ev, epsmin, A, x, b, lambda);
%     plot(cost(:));title('cost');
%     drawnow;
    if i >= 2
        if abs(cost(i,1)-cost(i-1,1)) / abs(cost(i,1)) <= 1e-5
            break;
        end
    end
    % If you have ground truth, compute SNR. If not comment these lines.
    SNR_iter = 20 * log10(norm(X0(:))/norm(Xr(:)-X0(:)));
    fprintf('iter: %d----> SNR: %6.4f \n', i, SNR_iter);
    SNR(i) = SNR_iter;
end
runtime = toc; %�㷨����ʱ��
end