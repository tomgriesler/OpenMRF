%% init pulseq
% basic SPI (spiral) readout
clear
seq_name = 'rosette_interleaved';

pulseq_scanner = 'Siemens_Prisma_3T_Technion';

% optional flags
flag_backup = 1; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

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

%% sequence params
% basic params
SPI.NR     = 299;             % [ ] number of repetitions
SPI.TR     = 0 *1e-3;       % [s] repetition time, 0 for minimization
SPI.Ndummy = 0;             % [ ] initial dummy loops

% slice excitation
SPI.exc_mode       = 'sinc';        % 'sinc' or 'sigpy_SLR'
SPI.exc_shape      = 'ex';          % only for sigpy: 'st' or 'ex' 
SPI.exc_time       = 2.0 *1e-3;     % [s] excitation time
SPI.exc_tbw        = 4;             % [ ] time bandwidth product
SPI.exc_fa_mode    = 'equal';       % 'equal',  'ramped',  'import'  
SPI.exc_flipangle  = 90 *pi/180;    % [rad] const FA or start of FA ramp -> set [] for auto mode

% gradient spoiling
SPI.spoil_nTwist   = 8;          % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 3.5 *1e-3;  % [s] time for spoiler and rewinder gradients

% rf spoiling
SPI.spoil_rf_mode = 'quad';       % rf spoiling mode: 'lin' or 'quad'
SPI.spoil_rf_inc  = 117 *pi/180;  % rf spoiling increment [rad]

%% spiral geoemtry mode
SPI.geo.traj_mode       = 'import';         % 'standard' or 'vds'
SPI.geo.interleave_mode = 'RandomGoldenAngle';    % 'Equal2Pi', 'Random2Pi', 'GoldenAngle', 'RandomGoldenAngle', 'RoundGoldenAngle'

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

SPI.Trec = 2;
SPI.Nrep = 1;

%% fat saturation and reset
FAT.mode = 'off';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% noise pre-scans
SPI.Nnoise = 100;
SPI_add_prescans();

%% create sequence

% saturtion or crusher
SAT_add();   

% recovery time
seq.addBlock(mr.makeDelay(SPI.Trec));

for loop_rep = 1 : SPI.Nrep
for loop_NR = 1-SPI.Ndummy : SPI.NR

    % saturtion or crusher
    SAT_add();   
    
    % recovery time
    seq.addBlock(mr.makeDelay(SPI.Trec));

    % fat saturation
    FAT_add();

    % spiral imaging
    SPI_add();

end
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();