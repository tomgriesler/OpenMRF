function [SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system, flag_plot)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ----------------- readout: spiral (SPI) -----------------
% ---------------------------------------------------------

% -------------------------------------------------------------------
% calculate rf pulses and slice selection gradient
% calculate slice rephaser (for 2D) or partition gradients (for 3D)
% calculate spiral k-space trajectory and gradients
% calculate slice spoiler or rewinder gradients
% set rf spoiling
% adapt timings for TRs and TEs
% -------------------------------------------------------------------

% optional import parameters for MRF use via SPI.mrf_import:
% phi, dphi, phi_id -> rotation angles of spiral interleaves
% FAs -> flip agnles
% TRs -> repetition times
% TEs -> echo times
% rf_phase -> rf spoiling increments
% gxy -> gradient trajecotry
% kxy -> k-space trajecotry

%% calculate rf pulses and slice selection gradient

% calculate or import flip angles
switch SPI.exc_fa_mode
    case 'equal'
        SPI.exc_flipangle = SPI.exc_flipangle * ones(SPI.NR,1);
    case 'ramped'
        SPI.exc_flipangle = SPI_calc_rf_ramp(SPI.NR, SPI.exc_flipangle);
    case 'import'
        SPI.exc_flipangle = SPI.mrf_import.FAs;
end

% default slice grad limit
if ~isfield(SPI, 'lim_gz_slew')
    SPI.lim_gz_slew = 1.0;
end
temp_sys = system;
temp_sys.maxSlew = temp_sys.maxSlew * SPI.lim_gz_slew;

% create rf and gz objects
for j=1:SPI.NR
    if strcmp(SPI.exc_mode, 'sinc')
        [SPI.rf(j,1), SPI.gz] = mr.makeSincPulse( SPI.exc_flipangle(j), ...
                                                  temp_sys, ...
                                                 'Duration', SPI.exc_time,...
                                                 'SliceThickness', FOV.dz, ...
                                                 'timeBwProduct', SPI.exc_tbw, ...
                                                 'apodization', 0.5, ...
                                                 'PhaseOffset', 0, ...
											     'use', 'excitation' );
    elseif strcmp(SPI.exc_mode, 'sigpy_SLR')
        [SPI.rf(j,1), SPI.gz] = SIGPY_SLR( SPI.exc_flipangle(j), ...
                                           SPI.exc_time, ...
                                           0, ...
                                           SPI.exc_tbw, ...
                                           SPI.exc_shape, ...
                                           'ls', 0.01, 0.01, 0 , ...
                                           FOV.dz, ...
										   temp_sys);
    end
    SPI.rf(j).freqOffset = FOV.z_offset * SPI.gz.amplitude;
end

clear temp_sys;

%% calculate slice rephaser (for 2D) or partition gradients (for 3D)

% check 2D or 3D mode
if FOV.Nz > 1
    SPI.mode_3D = 'on';
else
    SPI.mode_3D = 'off';
end

% default limit
if ~isfield(SPI, 'lim_reph_grad')
    SPI.lim_reph_grad = 0.9;
end
if ~isfield(SPI, 'lim_reph_slew')
    SPI.lim_reph_slew = 0.9;
end

% calculate slice rephaser for 2D
if strcmp(SPI.mode_3D, 'off')
    SPI.gz_reph(1)    = mr.makeTrapezoid('z', 'Area', -SPI.gz.area/2, 'maxGrad', system.maxGrad*SPI.lim_reph_grad, 'maxSlew', system.maxSlew*SPI.lim_reph_slew, 'system', system);
    SPI.reph_duration = mr.calcDuration(SPI.gz_reph(1));
end

% calculate slice rephaser for 3D
if strcmp(SPI.mode_3D, 'on')
    SPI.deltakz = 1 / FOV.fov_z;
    SPI.kz_area = (-FOV.Nz/2 : 1 : FOV.Nz/2-1)' * SPI.deltakz;
    for j = 1:FOV.Nz
        SPI.gz_reph(j,1) = mr.makeTrapezoid('z', 'Area', -SPI.gz.area/2 + SPI.kz_area(j), 'maxGrad', system.maxGrad*SPI.lim_reph_grad, 'maxSlew', system.maxSlew*SPI.lim_reph_slew, 'system', system);
        temp_dur(j)      = mr.calcDuration(SPI.gz_reph(j));
    end
    SPI.reph_duration = max(temp_dur); clear temp_dur;
    for j = 1:FOV.Nz
        SPI.gz_reph(j,1) = mr.makeTrapezoid('z', 'Area', -SPI.gz.area/2 + SPI.kz_area(j), 'Duration', SPI.reph_duration, 'maxGrad', system.maxGrad*SPI.lim_reph_grad, 'maxSlew', system.maxSlew*SPI.lim_reph_slew, 'system', system);
    end    
end
clear j;

%% calculate spiral k-space trajectory and gradients

% calculate basic k-space geometry
switch SPI.geo.traj_mode
    case 'standard'
        SPI.geo.traj = SPI_traj_standard(SPI.geo.kmax, SPI.geo.adcTime, SPI.geo.N, SPI.geo.a, SPI.geo.b, system.gradRasterTime, SPI.geo.adcNSamples);
    case 'vds'
        SPI.geo.traj = SPI_traj_vds(SPI, FOV, system, 0);
    case 'import'
        SPI.geo.traj = SPI_traj_import(SPI, system);
end

% calculate ADC
[SPI.adc, SPI.adcTime, SPI.adcNSamples, SPI.adcBW, SPI.adcDwell] = SPI_calc_adc(SPI, system);

% adc padding: start the gx gy gradients after adc dead time
SPI.adcTPad = ceil(  5   * system.adcDeadTime / system.gradRasterTime ) * system.gradRasterTime; % no gx and gy during start of adc
SPI.adcNPad = round( 2.5 * system.adcDeadTime / SPI.adc.dwell );                                 % delete samples at start of adc

% calculate or import rotation angles
switch SPI.geo.interleave_mode
    case 'Equal2Pi'
        SPI.phi_id = 1 : SPI.NR;
        SPI.dphi   = 2*pi / SPI.NR;
        SPI.phi    = SPI.phi_id * SPI.dphi;
    case 'Random2Pi'
        [~, SPI.phi_id] = sort(rand(1,SPI.NR));
        SPI.dphi        = 2*pi / SPI.NR;
        SPI.phi         = SPI.phi_id * SPI.dphi;
    case 'GoldenAngle'
        SPI.phi_id = 1 : SPI.NR;
        SPI.dphi   = 4*pi / (1+sqrt(5));
        SPI.phi    = SPI.phi_id * SPI.dphi;
    case 'RandomGoldenAngle'
        [~, SPI.phi_id] = sort(rand(1,SPI.NR));
        SPI.dphi        = 4*pi / (1+sqrt(5));
        SPI.phi         = SPI.phi_id * SPI.dphi;
    case 'RoundGoldenAngle'
        SPI.dphi   = 4*pi / (1+sqrt(5));
        SPI.phi    = wrapTo2Pi(SPI.dphi * (1:SPI.NR));
        SPI.dphi   = 2*pi / SPI.Nunique;
        SPI.phi    = round(SPI.phi/SPI.dphi)' * SPI.dphi;
        SPI.phi_id = round(SPI.phi/SPI.dphi);
        SPI.phi_id(SPI.phi_id==0) = SPI.Nunique;
        SPI.phi = SPI.phi_id * SPI.dphi;
    case 'import'
        SPI.phi    = SPI.mrf_import.phi;
        SPI.phi_id = SPI.mrf_import.phi_id;
end
if size(SPI.phi_id,1)<size(SPI.phi_id,2)
    SPI.phi_id = SPI.phi_id';
end
if size(SPI.phi,1)<size(SPI.phi,2)
    SPI.phi = SPI.phi';
end

% start with 1st golden angle
SPI.phi = wrapTo2Pi(SPI.phi - SPI.phi(1) + 2*pi/(1+sqrt(5))); 

% build unique phi array and check consistency
SPI.phi_unique = unique([SPI.phi_id, SPI.phi], 'rows');
SPI.phi_unique(:,1) = [];
if isfield(SPI, 'Nunique')
    if SPI.Nunique ~= numel(SPI.phi_unique)
        error('inconsistent number of unique spirals');
    end
else
    SPI.Nunique = numel(SPI.phi_unique);
end

% calculate spiral gradients including rewinders and spoilers
for j = 1:SPI.Nunique
    [SPI.gx(j,1), SPI.gy(j,:), SPI.sx(j,:), SPI.sy(j,:)] = SPI_calc_gx_gy(SPI, SPI.phi_unique(j), system);
end

%% calculate slice spoiler or rewinder gradients

% check balanced or unbalanced mode
if abs(SPI.spoil_nTwist) > 0
    SPI.spoil_gz_mode = 'unbalanced';
else
    SPI.spoil_gz_mode = 'balanced';
end

% default limits
if ~isfield(SPI, 'lim_spoil_grad')
    SPI.lim_spoil_grad = 1/sqrt(3);
end
if ~isfield(SPI, 'lim_spoil_slew')
    SPI.lim_spoil_slew = 1/sqrt(3);
end

% final kz value at end of TR
if strcmp(SPI.spoil_gz_mode, 'balanced')
    SPI.spoil_nTwist = 0;
end
kz_final = SPI.spoil_nTwist / FOV.dz;

% calculate spoiler/rewinder for 2D
if strcmp(SPI.mode_3D, 'off')
    SPI.gz_spoil(1) = mr.makeTrapezoid('z', 'Area', kz_final - SPI.gz.area/2, 'Duration', SPI.spoil_duration, 'maxGrad', system.maxGrad*SPI.lim_spoil_grad, 'maxSlew', system.maxSlew*SPI.lim_spoil_slew, 'system', system);  
end

% calculate spoiler/rewinder for 3D
if strcmp(SPI.mode_3D, 'on')
    for j = 1:FOV.Nz
        SPI.gz_spoil(j,1) = mr.makeTrapezoid('z', 'Area', kz_final - SPI.gz.area/2 - SPI.kz_area(j), 'Duration', SPI.spoil_duration, 'maxGrad', system.maxGrad*SPI.lim_spoil_grad, 'maxSlew', system.maxSlew*SPI.lim_spoil_slew, 'system', system);
    end    
end

% add spoiler delay
if ~isfield(SPI, 'spoil_gz_timing')
    SPI.spoil_gz_timing = 'fast';
end
if strcmp(SPI.spoil_gz_timing, 'fast')
    SPI.spoil_delay = mr.calcDuration(SPI.adc) + 2*system.gradRasterTime;
    SPI.spoil_delay = ceil(SPI.spoil_delay/system.gradRasterTime) * system.gradRasterTime;
end
if strcmp(SPI.spoil_gz_timing, 'slow')
    SPI.spoil_delay = mr.calcDuration(SPI.gx, SPI.gy, SPI.adc) + 2*system.gradRasterTime;
    SPI.spoil_delay = ceil(SPI.spoil_delay/system.gradRasterTime) * system.gradRasterTime;
end
for j = 1:FOV.Nz
    SPI.gz_spoil(j).delay = SPI.spoil_delay;
end 

clear kz_final j;

%% set rf spoiling
switch SPI.spoil_rf_mode
    case 'quad'
        SPI.spoil_rf_pow = 2;
    case 'lin'
        SPI.spoil_rf_pow = 1;
    case 'import'
        SPI.spoil_rf_pow = []; 
    otherwise
        error('unknown rf spoiling mode!')
end

%% adapt timings for TRs and TEs

% create or import TR list
if ~isfield(SPI, 'TR')
    SPI.TR = zeros(SPI.NR,1);
end
if isfield(SPI, 'mrf_import')
    if isfield(SPI.mrf_import, 'TRs')
        if numel(SPI.mrf_import.TRs) == SPI.NR
            SPI.TR = SPI.mrf_import.TRs;
        else
            error('incorrect number of TRs');
        end
    end
else
    if numel(SPI.TR) < SPI.NR
        SPI.TR = SPI.TR * ones(SPI.NR,1);
    end
end
if size(SPI.TR,1)<size(SPI.TR,2)
    SPI.TR = SPI.TR';
end

% create or import TE list
if ~isfield(SPI, 'TE')
    SPI.TE = zeros(SPI.NR,1);
end
if isfield(SPI, 'mrf_import')
    if isfield(SPI.mrf_import, 'TEs')
        if numel(SPI.mrf_import.TEs) == SPI.NR
            SPI.TE = SPI.mrf_import.TEs;
        else
            error('incorrect number of TEs');
        end
    end
else
    if numel(SPI.TE) < SPI.NR
        SPI.TE = SPI.TE * ones(SPI.NR,1);
    end
end
if size(SPI.TE,1)<size(SPI.TE,2)
    SPI.TE = SPI.TE';
end

% calculate TE filling delays
temp_te_min  = mr.calcDuration(SPI.gz)/2 + mr.calcDuration(SPI.gz_reph(1)) + system.gradRasterTime;
temp_te_fill = SPI.TE - temp_te_min;
if sum(temp_te_fill<0)>0
    warning('negative TE filling delay! TE corrected!');
end
for j = 1:SPI.NR
    if temp_te_fill(j) >= system.gradRasterTime
        SPI.TE_delay(j) = mr.makeDelay( round(temp_te_fill(j)/system.gradRasterTime) * system.gradRasterTime );
    else
        SPI.TE_delay(j,1) = mr.makeDelay( system.gradRasterTime );
        SPI.TE(j)         = temp_te_min;
    end
end
clear temp_te_min temp_te_fill j;

% calculate TR filling delay
for j = 1:SPI.NR
    temp_tr_min(j,1) = mr.calcDuration(SPI.rf, SPI.gz) + mr.calcDuration(SPI.gz_reph(1)) + mr.calcDuration(SPI.TE_delay(j)) + mr.calcDuration(SPI.gx, SPI.gy, SPI.gz_spoil, SPI.adc) + system.gradRasterTime;
end
temp_tr_fill = SPI.TR - temp_tr_min;
if sum(temp_tr_fill<0)>0
    warning('negative TR filling delay! TR corrected!');
end
for j = 1:SPI.NR
    if temp_tr_fill(j) >= system.gradRasterTime
        SPI.TR_delay(j,1) = mr.makeDelay( round(temp_tr_fill(j)/system.gradRasterTime) * system.gradRasterTime );
    else
        SPI.TR_delay(j,1) = mr.makeDelay( system.gradRasterTime );
        SPI.TR(j)         = temp_tr_min(j);
    end
end
clear temp_tr_min temp_tr_fill j;

if ~isfield(SPI, 'Ndummy')
    SPI.Ndummy = 0;
end

%% export kspace trajectory for reconstruction
seq = mr.Sequence(system);
if nargin<4
    flag_plot = 0;
end
if isfield(SPI, 'kz_area')
    loop_kz = find(SPI.kz_area==0); % only center partition; reduce array size and computation time
else
    loop_kz = 1;
end
for loop_NR = 1:SPI.Nunique
    seq.addBlock(SPI.rf(1), SPI.gz);
    seq.addBlock(SPI.gz_reph(loop_kz));
    seq.addBlock(SPI.gx(loop_NR), SPI.gy(loop_NR), SPI.gz_spoil(loop_kz), SPI.adc);  % adc loop
end
[ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, flag_plot);
tempx = ktraj_adc(1,:);
tempy = ktraj_adc(2,:);
if size(tempx,2) > SPI.Nunique*SPI.adcNSamples
    tempx(SPI.Nunique*SPI.adcNSamples+1:end) = [];
    tempy(SPI.Nunique*SPI.adcNSamples+1:end) = [];
end
tempx = reshape(tempx', [numel(tempx)/SPI.Nunique,SPI.Nunique])';
tempy = reshape(tempy', [numel(tempy)/SPI.Nunique,SPI.Nunique])';
ktraj_reco(1,:,:) = tempx(:,:);
ktraj_reco(2,:,:) = tempy(:,:);
ktraj_adc  = []; % reduce backup size
ktraj_full = []; % reduce backup size
clear tempx tempy loop_NR t_excitation t_refocusing;

end
