% reconstruction of Siemens inversion-recovery (IR) spin-echo (SE) or turbo-spin-echo (TSE)
% use for T1 mapping of the NIST phantom
% mgram; V1; 28.11.2024

%% start reco
clear

% read rawdata from .dat
study_path = 'Q:/data/Pulseq/Rawdata/tomgr/CimaX_Freiburg/250823_Open_MRF__NIST_T1_T2_T1p';
study_name = { 'meas_MID00373_FID31834_t1_se_cor_TI_35ms.dat',
               'meas_MID00374_FID31835_t1_se_cor_TI_60ms.dat',
               'meas_MID00375_FID31836_t1_se_cor_TI_95ms.dat',
               'meas_MID00376_FID31837_t1_se_cor_TI_160ms.dat',
               'meas_MID00377_FID31838_t1_se_cor_TI_260ms.dat',
               'meas_MID00378_FID31839_t1_se_cor_TI_430ms.dat',
               'meas_MID00379_FID31840_t1_se_cor_TI_710ms.dat',
               'meas_MID00380_FID31841_t1_se_cor_TI_1180ms.dat',
               'meas_MID00381_FID31842_t1_se_cor_TI_1940ms.dat',
               'meas_MID00382_FID31843_t1_se_cor_TI_3200ms.dat'};

for j=1:numel(study_name)
    twix_obj = mapVBVD(fullfile(study_path, study_name{j}), 'ignoreSeg', 'removeOS');
    twix_obj = twix_obj{end};
    rawdata(j,:,:,:) = permute(twix_obj.image(),[2,3,1]);
end

% ifft
Images_coils = kspace2image(rawdata);

% cmaps: use the last image as the reference
[cmaps, ~, ~] = mg_espirit_cmaps(squeeze(Images_coils(end,:,:,:)), [], [], [], []);

% coil combined images
[NInv, NCoils, Ny, Nx] = size(Images_coils);
Images = zeros(NInv, Ny, Nx);
for j=1:NInv
    Images(j,:,:,:) = squeeze(sum(squeeze(Images_coils(j,:,:,:)) .* conj(cmaps)));
end

% zero interpolation filling
zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;
Images = mg_zero_filling(Images, zero_params);

% fit mask
[mask, mask3D] = mg_get_mask_fit(squeeze(mean(abs(Images))), 'holes', NInv);

%% T1 mapping

% rotate images to real axis: use the last image as the reference
Images = real(Images .* exp(-1i*repmat(angle(Images(end,:,:)),[NInv,1,1])));

TI = [35 60 95 160 260 430 710 1180 1940 3200] *1e-3;
[T1_Map, M0_Map, Eff_Map, R2_Map] = mg_map_T1(real(Images), TI, mask);

%% vis results
t1lims = [0 2500] *1e-3;
t1cmp  = get_cmp('T1', 1000, 1);

figure()
ax1 = subplot(1,3,1);
imagesc(T1_Map, t1lims); axis image; axis off; colormap(ax1, t1cmp); colorbar; title('T1 map [s]');
ax2 = subplot(1,3,2);
imagesc(R2_Map, [0.8 1]); axis image; axis off; colormap(ax2, turbo(1000)); colorbar; title('R2 map');
ax3 = subplot(1,3,3);
imagesc(Eff_Map, [0.5, 1]); axis image; axis off; colormap(ax3, turbo(1000)); colorbar; title('Inv Eff Map');