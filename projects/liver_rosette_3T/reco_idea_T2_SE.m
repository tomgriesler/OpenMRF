% reconstruction of Siemens spin-echo (SE)
% use for T2 mapping of the NIST phantom
% mgram; V1; 28.11.2024

%% start reco
clear

% read rawdata from .dat
study_path = '';
study_name = { 'meas_MID00157_FID10912_se_te10.dat',
               'meas_MID00158_FID10913_se_te15.dat',
               'meas_MID00159_FID10914_se_te20.dat',
               'meas_MID00160_FID10915_se_te30.dat',
               'meas_MID00161_FID10916_se_te40.dat',
               'meas_MID00162_FID10917_se_te60.dat',
               'meas_MID00163_FID10918_se_te85.dat',
               'meas_MID00164_FID10919_se_te120.dat',
               'meas_MID00165_FID10920_se_te175.dat',
               'meas_MID00166_FID10921_se_te250.dat',
               'meas_MID00167_FID10922_se_te400.dat'};

for j=1:numel(study_name)
    twix_obj = mapVBVD(fullfile(study_path, study_name{j}), 'ignoreSeg', 'removeOS');
    twix_obj = twix_obj{end};
    twix_obj = squeeze(twix_obj.image());
    rawdata(j,:,:,:) = permute(twix_obj,[2,3,1]);
end

% ifft
Images_coils = kspace2image(rawdata);

% cmaps: use the first image as the reference
[cmaps, ~, ~] = mg_espirit_cmaps(squeeze(Images_coils(1,:,:,:)), [], [], [], []);

% coil combined images
[NInv, NCoils, Ny, Nx] = size(Images_coils);
Images = zeros(NInv, Ny, Nx);
for j=1:NInv
    Images(j,:,:,:) = squeeze(sum(squeeze(Images_coils(j,:,:,:)) .* conj(cmaps)));
end

% zero interpolation filling
zero_params.onoff  = 1;
zero_params.radius = 0.5;
zero_params.factor = 2.0;
Images = mg_zero_filling(Images, zero_params);

% fit mask
[mask, mask3D] = mg_get_mask_fit(squeeze(mean(abs(Images))), 'holes', NInv);

%% T2 mapping

% rotate images to real axis: use the first image as the reference
Images = real(Images .* exp(-1i*repmat(angle(Images(1,:,:)),[NInv,1,1])));

TE = [10 15 20 30 40 60 85 120 175 250 400] *1e-3;
[ T2_Map, ~, R2_Map ] = mg_map_T12p( Images, TE, mask );

%% vis results
t2lims = [0 1000] *1e-3;
t2cmp  = get_cmp('T2', 1000, 1);

figure()
ax1 = subplot(1,2,1);
imagesc(T2_Map, t2lims); axis image; axis off; colormap(ax1, t2cmp); colorbar; title('T2 map [s]');
ax2 = subplot(1,2,2);
imagesc(R2_Map, [0.8 1]); axis image; axis off; colormap(ax2, turbo(1000)); colorbar; title('R2 map');

%%
res.T2_Map = T2_Map;
res.R2_Map = R2_Map;
res.mask = mask;
res.Images = Images;
res.TE = TE;

save_study_results('t2_se_cartesian', res, study_path);