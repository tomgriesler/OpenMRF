function phantom = mg_numerical_phantom(Nxy, Ncoils, coil_params, geo)

%% define phantom basics 

% Nxy:    grid size
% Ncoils: number of simulated coils

% default coil params
if isempty(coil_params)
    coil_params.Nb_coils   = Ncoils;
    coil_params.res        = Nxy * [1, 1];
    coil_params.type       = 'biot';
    coil_params.param.FS   = 0.28; % FOV width
    coil_params.param.D    = 0.25; % Distance center->coil
    coil_params.param.R    = 0.2;  % radius of the coil
    coil_params.sens_model = 'sinusoidal';
end

% geoemtry similar to the Caliber MRI phantom
if isempty(geo)
    n_inner = 5;    % number of circular phantom inserts (inner)
    n_outer = 9;    % number of circular phantom inserts (outer)
    r_phant = 0.9;  % radius of main phantom
    r_outer = 0.65; % radius of circle along which outer phantoms lay
    r_inner = 0.3;  % radius of circle along which inner phantoms lay
    r_insrt = 0.1;  % radius of each circular insert
else
    n_inner = geo.n_inner;
    n_outer = geo.n_outer;
    r_phant = geo.r_phant;
    r_outer = geo.r_outer;
    r_inner = geo.r_inner;
    r_insrt = geo.r_insrt;
end

% total number of circular phantom inserts
n_phant = n_inner + n_outer;

%% simulate coils
% https://bigwww.epfl.ch/algorithms/mriphantom/
% https://bigwww.epfl.ch/algorithms/mriphantom/MRIPhantomv0-8.zip
% Guerquin-Kern M, Lejeune L, Pruessmann KP, Unser M.
% Realistic analytical phantoms for parallel magnetic resonance imaging.
% IEEE Trans Med Imaging. 2012 Mar;31(3):626-36.
% doi: 10.1109/TMI.2011.2174158

% init basic circular phantom
circ.FOV = [1, 1];
circ.region{1}.type   = 'ellipse';
circ.region{1}.center = [0, 0];
circ.region{1}.angle  = 0;
circ.region{1}.weight = 1;
circ.region{1}.width  = r_phant * [1, 1];
res                   = Nxy * [1, 1];
mask2D                = imfill(RasterizePhantom(circ, res, 1)==1, 'holes');

% add NIST-like edges
% [x, y] = meshgrid(1:Nxy, 1:Nxy);
% phi = linspace(0, 2*pi, 9) + pi/7;
% phi(end) = [];
% for j=1:8
%     xc = round(Nxy/2 + r_phant*Nxy/2*1.08 * cos(phi(j)));
%     yc = round(Nxy/2 + r_phant*Nxy/2*1.08 * sin(phi(j)));
%     mask2D(((x - xc).^2 + (y - yc).^2) <= (0.08*Nxy)^2) = false;
% end
% clear x y xc yc phi;

% parameter of the model: polynomial degree or bandwith
if strcmp(coil_params.sens_model, 'polynomial')
    param = 8;
end
if strcmp(coil_params.sens_model, 'sinusoidal')
    param = 6;
end
coil_params = simulate_sensitivities(coil_params);

% coil simultion
cmaps    = zeros(Ncoils, Nxy, Nxy);
s        = cell(1, Ncoils);
sens_num = cell(1, Ncoils);
residue  = cell(1, Ncoils);
for j = 1 : Ncoils
    sens             = coil_params.sensitivity(:,:,j);
    sens             = sens / max(sens(mask2D));
    s{j}             = SensFitting(sens, coil_params.sens_model, param, mask2D);
    [~, sens_num{j}] = RasterizePhantom(circ, res, s{j});
    residue{j}       = sens_num{j}(mask2D) - sens(mask2D);
    cmaps(j,:,:)     = sens_num{j} .* mask2D;
    clear sens;
end

cmaps = cmaps / max(abs(cmaps(:)));

% figure('Name', 'abs(cmaps)', 'NumberTitle', 'off');
% mycmp = jet(1000);
% mycmp(1,:) = 0;
% for j=1:Ncoils
%     subplot(1, Ncoils, j)
%     imagesc(abs(squeeze(cmaps(j,:,:))) + 10*(mask2D-1), [0, 1]);   axis image; axis off; colormap(mycmp);
%     sgtitle('abs() of cmaps')
% end

% figure('Name', 'angle(cmaps)', 'NumberTitle', 'off');
% mycmp = plasma(1000);
% mycmp(1,:) = 0;
% for j=1:Ncoils
%     subplot(1, Ncoils, j)
%     imagesc(angle(squeeze(cmaps(j,:,:))) + 10*(mask2D-1), max(abs(angle(cmaps(:))))*[-1 1]); axis image; axis off; colormap(mycmp);
%     sgtitle('angle() of cmaps')
% end

%% create circular phantom inserts -> ind_map 

% angles of circular inserts
phi_inner = linspace(0, 2*pi, n_inner+1); phi_inner(end) = []; % inner circle
phi_outer = linspace(0, 2*pi, n_outer+1); phi_outer(end) = []; % outer circle
phi_inner = phi_inner + pi/3; % twist inner circle w.r.t. outer circle
phi_outer = phi_outer + pi/5; % twist outer circle w.r.t. inner circle

% inner circles
ind_map  = zeros(Nxy,Nxy);
xC_inner = r_inner * cos(phi_inner); xC_inner = xC_inner*(Nxy/2) + Nxy/2;
yC_inner = r_inner * sin(phi_inner); yC_inner = -yC_inner*(Nxy/2) + Nxy/2;
rC_inner = r_insrt * ones(length(xC_inner),1)'; rC_inner = rC_inner * (Nxy/2);
circles_inner = transpose([xC_inner; yC_inner; rC_inner]);
for j=1:length(circles_inner)
    for y=1:Nxy
    for x=1:Nxy
        if ( (x-circles_inner(j,1))^2 + (y-circles_inner(j,2))^2 ...
                < circles_inner(j,3)^2)
            ind_map(y,x) = j+1;
        end
    end
    end
end
clear j x y

% outer circles 
xC_outer = r_outer * cos(phi_outer); xC_outer =  xC_outer*(Nxy/2) + Nxy/2;
yC_outer = r_outer * sin(phi_outer); yC_outer = -yC_outer*(Nxy/2) + Nxy/2;
rC_outer = r_insrt * ones(length(xC_outer),1)'; rC_outer = rC_outer * (Nxy/2);
circles_outer = transpose([xC_outer; yC_outer; rC_outer]);
for j=1:length(circles_outer)
    for y=1:Nxy
    for x=1:Nxy
        if ( (x-circles_outer(j,1))^2 + (y-circles_outer(j,2))^2 ...
                < circles_outer(j,3)^2)
            ind_map(y,x) = j+n_inner+1;
        end
    end
    end
end
clear j x y

% background -> index 1
ind_map((mask2D==1) & (ind_map==0)) = 1;

figure('Name', 'Index map', 'NumberTitle', 'off');
imagesc(ind_map); axis image; axis off; colormap(parula);
title('index map')

%% randomize receiver phase offsets
rng(1, "twister");
cmaps = cmaps.*exp(1i*repmat(rand(size(cmaps,1),1)*2*pi, 1, size(cmaps,2), size(cmaps,3)));

%% output phantom as struct
phantom.Nxy         = Nxy;
phantom.Ncoils      = Ncoils;
phantom.n_phant     = n_phant;
phantom.geo.n_inner = n_inner;
phantom.geo.n_outer = n_outer;
phantom.geo.r_phant = r_phant;
phantom.geo.r_inner = r_inner;
phantom.geo.r_outer = r_outer;
phantom.geo.r_insrt = r_insrt;
phantom.coil_params = coil_params;
phantom.cmaps       = cmaps;
phantom.ind_map     = ind_map;
phantom.mask2D      = mask2D;
phantom.mask1D      = mask2D(:)==1;

end
