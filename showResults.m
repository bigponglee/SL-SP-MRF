function [] = showResults(T1, T2, PD, k_space_undersampling_ratio)
%Compare results between methods: T1, T2 and PD maps
% close all;

%% T1
T1_value = 2500;
T1_value_error = 500;
T2_value = 500;
T2_value_error = 200;
PD_value = 116.8877;
PD_value_error = 23.3775;
figure('Name', 'T1 T2 PD weighted images', 'NumberTitle', 'off');
colormap(hot);
subplot(3, 3, 1);
imagesc(T1.full, [0, T1_value]);
title('T1 Gold standard', 'FontSize', 12.5);
axis off;
axis image;
colorbar;
subplot(3, 3, 2);
imagesc(T1.old_mrf, [0, T1_value]);
% title(['Fully sampled MRF'],'FontSize',12.5);axis off;
title(['MRF ', num2str(k_space_undersampling_ratio*100), '% sampled'], 'FontSize', 12.5);
axis off;
axis image;
colorbar;
subplot(3, 3, 3);
imagesc(T1.giraf, [0, T1_value]);
% set(gca,'fontsize',12.5);
title(['GIRAF ', num2str(k_space_undersampling_ratio*100), '% sampled'], 'FontSize', 12.5);
axis off;
axis image;
colorbar;

%% T2
subplot(3, 3, 4);
imagesc(T2.full, [0, T2_value]);
title('T2 Gold standard', 'FontSize', 12.5);
axis off;
axis image;
colorbar;
subplot(3, 3, 5);
imagesc(T2.old_mrf, [0, T2_value]);
axis off;
axis image;
colorbar;
subplot(3, 3, 6);
imagesc(T2.giraf, [0, T2_value]);
axis off;
axis image;
colorbar;

%% PD
subplot(3, 3, 7);
imagesc(PD.full, [0, PD_value]);
title('PD Gold standard', 'FontSize', 12.5);
axis off;
axis image;
colorbar;
subplot(3, 3, 8);
imagesc(PD.old_mrf, [0, PD_value]);
axis off;
axis image;
colorbar;
subplot(3, 3, 9);
imagesc(PD.giraf, [0, PD_value]);
axis off;
axis image;
colorbar;

%% Error maps
figure('Name', 'Error maps', 'NumberTitle', 'off');
colormap(hot);

%% T1
mse_old_mrf = goodnessOfFit(T1.old_mrf(:), T1.full(:), 'NMSE');
subplot(3, 2, 1);
imagesc(abs(T1.old_mrf-T1.full), [0, T1_value_error]);
axis off;
axis image;
colorbar;
% title('Fully Sampled MRF','FontSize',12.5);
title('MRF', 'FontSize', 20);
text(-20, 55, 'T1', 'FontSize', 15);
text(20, 140, ['NMSE=', num2str(mse_old_mrf)], 'FontSize', 12.5);

mse_cs = goodnessOfFit(T1.giraf(:), T1.full(:), 'NMSE');
subplot(3, 2, 2);
imagesc(abs(T1.giraf-T1.full), [0, T1_value_error]);
axis off;axis image;
% colorbar;set(gca,'fontsize',12.5);
title('GIRAF', 'FontSize', 20);
text(20, 140, ['NMSE=', num2str(mse_cs)], 'FontSize', 12.5);
colorbar;

%% T2
mse_old_mrf = goodnessOfFit(T2.old_mrf(:), T2.full(:), 'NMSE');
subplot(3, 2, 3);
imagesc(abs(T2.old_mrf-T2.full), [0, T2_value_error]);
axis off;
axis image;
text(-20, 55, 'T2', 'FontSize', 15);
text(20, 140, ['NMSE=', num2str(mse_old_mrf)], 'FontSize', 12.5);
colorbar;

mse_cs = goodnessOfFit(T2.giraf(:), T2.full(:), 'NMSE');
subplot(3, 2, 4);
imagesc(abs(T2.giraf-T2.full), [0, T2_value_error]);
axis off;axis image;
text(20, 140, ['NMSE=', num2str(mse_cs)], 'FontSize', 12.5);
colorbar;

%% PD
mse_old_mrf = goodnessOfFit(PD.old_mrf(:), PD.full(:), 'NMSE');
subplot(3, 2, 5);
imagesc(abs(PD.old_mrf-PD.full), [0, PD_value_error]);
axis off;
axis image;
text(-22, 55, 'PD', 'FontSize', 15);
text(20, 140, ['NMSE=', num2str(mse_old_mrf)], 'FontSize', 12.5);
colorbar;

mse_cs = goodnessOfFit(PD.giraf(:), PD.full(:), 'NMSE');
subplot(3, 2, 6);
imagesc(abs(PD.giraf-PD.full), [0, PD_value_error]);
axis off;axis image;
text(20, 140, ['NMSE=', num2str(mse_cs)], 'FontSize', 12.5);
colorbar;

end
