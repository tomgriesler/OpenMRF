%% init pulseq
% basis: SPI readout
% use for: cardiac fingerprinting
clear
seq_name = 'amrf_t1_t2_3T';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

% select scanner
pulseq_scanner = 'Siemens_Prisma_3T_Technion';

% select pns sim orientation
pns_orientation = 'coronal';

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nxy      = 256;         % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 400  *1e-3;  % [m] FOV geometry
FOV.dz       = 8   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% MRF contrast encoding

% encoding list with contrast preparations for different segments
% 'No_Prep'    ->  use for recovery of longitudinal magnetization
% 'Saturation' ->  use for T1 encoding
% 'Inversion'  ->  use for T1 encoding
% 'T2'         ->  use for T2 encoding
% 'SL'         ->  use for T1p encoding
% 'ADIASL'     ->  use for adiabatic T1p encoding
% 'MLEV'       ->  use for T2 or T2p encoding

MRF.enc_list = {
'Inversion';
'No_Prep';
'T2';
'T2';
'T2';
'Inversion';
'No_Prep';
'T2';
'T2';
'T2';
'Inversion';
'No_Prep';
'T2';
'T2';
'T2';
'Inversion';
'No_Prep';
'T2';
'T2';
'T2';
};

MRF.n_segm = numel(MRF.enc_list);

%% MRF flipangle and repetition times: variable FA
MRF.nr  = 30;                   % numer of readouts per hear beat
MRF.NR  = MRF.n_segm * MRF.nr;  % total number of readouts
% MRF.FAs = deg2rad(15)*ones(MRF.NR, 1);
MRF.FAs = repmat(deg2rad([5:1.25:13.75 ones(1, 22)*15]), 1, MRF.n_segm);
MRF.TRs = 17 *1e-3 *ones(MRF.NR,1); 
MRF.TRs(MRF.nr:MRF.nr:end-1) = MRF.TRs(MRF.nr:MRF.nr:end-1)+0.4;
% MRF.TRs = zeros(MRF.NR, 1);

%% params: Spiral Readouts

% import params from MRF struct
SPI.nr             = MRF.nr;    % [ ] number of readouts per heart beat
SPI.NR             = MRF.NR;    % [ ] number of repetitions
SPI.mrf_import.TRs = MRF.TRs;   % [s] repetition times
SPI.mrf_import.FAs = MRF.FAs;   % [rad] flip angles

% slice excitation
SPI.exc_mode       = 'sinc';      % 'sinc' or 'sigpy_SLR'
SPI.exc_shape      = 'ex';        % only for sigpy: 'st' or 'ex' 
SPI.exc_time       = 0.8 *1e-3;   % [s] excitation time
SPI.exc_tbw        = 2;           % [ ] time bandwidth product
SPI.exc_fa_mode    = 'import';    % 'equal',  'ramped',  'import' 
SPI.lim_gz_slew    = 0.9;         % [ ] reduce stimulation during slice excitation
SPI.lim_reph_slew  = 0.9;         % [ ] reduce stimulation during slice rephaser

% gradient spoiling
SPI.spoil_nTwist   = 4;          % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 0.8 *1e-3;  % [s] time for spoiler and rewinder gradients !!!! 1.0ms
SPI.lim_spoil_slew = 0.9;        % [ ] reduce stimulation during gradient spoiling

% rf spoiling
SPI.spoil_rf_mode = 'lin';      % rf spoiling mode: 'lin' or 'quad'
SPI.spoil_rf_inc  = 0 *pi/180;  % rf spoiling increment [rad]

% spiral geometry mode
SPI.geo.interleave_mode = 'RoundGoldenAngle';
SPI.geo.traj_mode       = 'import';

% rosette geometry: import from poet
SPI.Nunique             = 299;                % number of unique projection angles
SPI.deltak              = 1/FOV.fov_xy;       % [1/m] kspace sampling
SPI.kmax                = SPI.deltak * FOV.Nxy/2;
% SPI.geo.path            = 'rosette_alt_21lobes_299arms_21gmx_400mm_256mat_240415.mat';
SPI.geo.path            = 'rosette_alt_17lobes_299arms_400mm_256mat_251108.mat';
SPI.geo.BW              = 400 *1e3;      % [Hz] bandwidth of spiral acquisition
SPI.geo.kmax            = SPI.kmax;      % determines resolution
SPI.geo.t_dwell         = 1/SPI.geo.BW;  % [s]  dwell time for spiral acquisition
SPI.spoil_gz_timing     = 'slow';

% calc rosette params
[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);
MRF.TRs = SPI.TR;

%% params: Inversion
INV.mode = 'on';
INV.rf_type = 'HYPSEC_inversion';
INV.tExc = 10 *1e-3;
INV.mu   = 4.9;
INV.beta = 700;
INV.inv_rec_time = [12 300 12 300] *1e-3;
INV = INV_init(INV, FOV, system);

%% params: T2 preparation
T2.rfc_dur    = 2 *1e-3;
T2.prep_times = [40 80 160 40 80 160 40 80 160 40 80 160] * 1e-3;
T2.exc_mode   = 'adiabatic_BIR4';
T2            = T2_init(T2, FOV, system);

%% check MRF encoding params
MRF_check_enc_list();

%% Trigger Mode
MRF.mode_trig = 'off';

for ii=1:MRF.n_segm
    MRF.delay_dynamic(ii, 1) = mr.makeDelay(1e-5);
end

MRF.delay_soft   = mr.makeDelay(1e-5); % minimize delays

%% ---------- add sequence objects ----------
MRF_add_segments();

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();