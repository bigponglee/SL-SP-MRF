function mat_of_contrasts = build_fully_sampled_contrasts(RFpulses, TR_times, TE, T1_image, T2_image, PD_image)
%This function simulates a fully sampled k-spaces for each echo with the MRF
%sequence and returns the fully sampled temporal contrasts
mat_of_contrasts = zeros(size(T1_image, 1), size(T1_image, 2), length(RFpulses), 'single');
locs = find((T2_image ~= 0) & (T1_image ~= 0));
T1_image = T1_image / 1000;
T2_image = T2_image / 1000;
TR_times = TR_times / 1000;
TE = TE / 1000;

for k = 1:length(locs)
    [i, j] = ind2sub(size(T1_image), locs(k));
    mat_of_contrasts(i, j, :) = (single((PD_image(i, j) * epg_fisp_mrf(RFpulses, TR_times, TE, T1_image(i, j), T2_image(i, j)))));
end
end
