close all;
clear;
addpath('simulations\');
addpath(genpath('methods\')); 
addpath('data\');
addpath('Bloch\');

%% 仿真数据
N = 128; %图像尺寸N*N
L = 500; %序列长度L
k_space_undersampling_ratio = 0.06; %采样率
gpuDevice(1); %GPU acc
load input_to_fisp_experiment %load the quantitative maps
[RFpulses, TR] = generate_RF_TR(L); %Load slowly changing RF and TR values
RFpulses = RFpulses * 1i; %to avoid complex values of X and D
TE = 2.94; %回波时间echo time（TE）
noisy = 0;

%% build dictionary using MRF FISP sequence, omit entries where T2>T1
if exist(['D_', num2str(L), '.mat'], 'file') ~= 0 %可以使用保存的字典、图像、LUT，减少时间
    load(['D_', num2str(L), '.mat'])
else
    disp('Building the dictionary...');
    [FISP_dictionary, LUT] = build_dictionary_fisp(L, RFpulses, TR, TE); %build the MRF FISP dictionary and its matching look up table
    D = single(FISP_dictionary);
    clear FISP_dictionary;
    LUT = LUT * 1000; %单位改为ms
    disp('The dictionary is ready, building the temporal images...');
    X = build_fully_sampled_contrasts(RFpulses, TR, TE, T1_128, T2_128, PD_128); %build the fully sampled temporal contrasts全采样信号
    disp('The images are ready');
    save(['.\data\D_', num2str(L), '.mat'], 'D', 'X', 'LUT')
end

%% 采样模板
load('radio_sampling_matrix_006.mat') %k_space_undersampling_ratio = 0.06
sampling_matrix = sampling_matrix(:, :, 1:L);

%% 噪声
kSpaceNoise = reshape([1, 1i]*0.5*randn(2, L*N^2), N, N, L); %高斯噪声

%% full sample
E_matched_filter_full = full(find_E_fast(X, D)); %results is one-sparse E
[T1_full, T2_full, PD_full] = build_maps_from_E(E_matched_filter_full, LUT);

%% GIRAF_SP
[X_estimated_giraf, Y_full, m_abs, time_GIRAF] = SL_SP(X, sampling_matrix, kSpaceNoise, D, noisy); %算法恢复矩阵
X_estimated_giraf = X_estimated_giraf * m_abs; %去归一化
E_matched_filter_giraf = full(find_E_fast(X_estimated_giraf, D)); %results is one-sparse E
[T1_giraf, T2_giraf, PD_giraf] = build_maps_from_E(E_matched_filter_giraf, LUT);
SNR_giraf = 20 * log10(norm(X(:))/norm(X_estimated_giraf(:)-X(:)));

%% build contrasts using conventional MRF
tic;
Y_full = Y_full * m_abs; %去归一化
Y = sampling_matrix .* Y_full; %under-sampling the noised data
Y = gather(reshape(Y, N, N*L));
Y_full = gather(reshape(Y_full, N, N*L));
X_estimated_old_mrf = reshape(fft_mats(Y, 2), N, N, L);
E_matched_filter_old_mrf = full(find_E_fast(X_estimated_old_mrf, D)); %results is one-sparse E
[T1_old_mrf, T2_old_mrf, PD_old_mrf] = build_maps_from_E(E_matched_filter_old_mrf, LUT);
SNR_old_mrf = 20 * log10(norm(X(:))/norm(X_estimated_old_mrf(:)-X(:)));
time_OldMRF = toc;

%% Show Results in the region of interest
T1.orig = T1_128;
T2.orig = T2_128;
PD.orig = PD_128;

T1.full = T1_full;
T2.full = T2_full;
PD.full = PD_full;

T1.old_mrf = T1_old_mrf;
T2.old_mrf = T2_old_mrf;
PD.old_mrf = PD_old_mrf;

T1.giraf = T1_giraf .* (T1.orig ~= 0);
T2.giraf = T2_giraf .* (T2.orig ~= 0);
PD.giraf = PD_giraf .* (PD.orig ~= 0);

showResults(T1, T2, PD, k_space_undersampling_ratio);