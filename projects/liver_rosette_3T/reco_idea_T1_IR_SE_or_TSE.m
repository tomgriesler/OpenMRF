% reconstruction of Siemens inversion-recovery (IR) spin-echo (SE) or turbo-spin-echo (TSE)
% use for T1 mapping of the NIST phantom
% mgram; V1; 28.11.2024

%% start reco
clear

% read rawdata from .dat
study_path = '';
study_name = { 'meas_MID00146_FID10901_ir_se_ti35.dat',
               'meas_MID00148_FID10903_ir_se_ti60.dat',
               'meas_MID00149_FID10904_ir_se_ti95.dat',
               'meas_MID00150_FID10905_ir_se_ti160.dat',
               'meas_MID00151_FID10906_ir_se_ti260.dat',
               'meas_MID00152_FID10907_ir_se_ti430.dat',
               'meas_MID00153_FID10908_ir_se_ti710.dat',
               'meas_MID00154_FID10909_ir_se_ti1180.dat',
               'meas_MID00155_FID10910_ir_se_ti1940.dat',
               'meas_MID00156_FID10911_ir_se_ti3200.dat'};

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

%%
res.T1_Map = T1_Map;
res.R2_Map = R2_Map;
res.Eff_Map = Eff_Map;
res.M0_Map = M0_Map;
res.mask = mask;
res.Images = Images;
res.TI = TI;

save_study_results('t1_ir_se_cartesian', res, study_path);