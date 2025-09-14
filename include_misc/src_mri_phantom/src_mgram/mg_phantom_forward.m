function DATA = mg_phantom_forward(Mxy, mask2D, cmaps, ktraj, phi_id)

% forward model for numerical phantom simulation
% Maximilian Gram, University of Würzburg, v1, 08.09.2025

% ----- input: -----
% Mxy:    NR x Nvoxel
% mask2D: Nxy x Nxy
% cmaps:  NCoils x Nxy x Nxy
% ktraj:  2 x NSpi x NRead
% phi_id: NR x 1

% ----- output: -----
% DATA: NCoils x NR x Nread

% get dimensions
[NR, ~]          = size(Mxy);
[NCoils, Nxy, ~] = size(cmaps);
NRead            = size(ktraj, 3);

% convert to single
Mxy    = single(Mxy);
cmaps  = permute(single(cmaps), [2,3,1]);
ktraj  = single(ktraj);
mask2D = mask2D==1;

% NUFFT operator
kx         = single(squeeze(ktraj(1,:,:)))';
ky         = single(squeeze(ktraj(2,:,:)))';
kmax       = max(sqrt(kx(:).^2 + ky(:).^2));
kx         = kx/kmax*Nxy/2;    
ky         = ky/kmax*Nxy/2;
FT         = NUFFT(kx/Nxy+1i*ky/Nxy, kx*0+1, [0 0], Nxy*[1 1]);
clear ktraj kx ky kmax fov Ni nufft_args G dcf_all;

% forward model:
% calc complex images for each TR -> add coils -> nufft -> apply sampling mask
DATA   = zeros(NR, NRead, NCoils, 'single');
idx_2D = reshape(1:(Nxy)^2, Nxy,Nxy);
idx_2D = idx_2D(mask2D);
for j = 1:NR
    temp_mxy    = Mxy(j,:);
    temp_image  = zeros(Nxy, Nxy);
    temp_image(idx_2D) = temp_mxy;
    temp_image  = repmat(temp_image, 1, 1, NCoils) .* cmaps;
    temp_signal = FT * temp_image;
    temp_signal = squeeze(temp_signal(:,phi_id(j),:));
    DATA(j,:,:) = temp_signal;
end
DATA = permute(DATA, [3, 1, 2]);

end

