%% init pulseq
% basis: SPI readout
% use for: T2 mapping
clear
seq_name = 't2_mapping_mlev_8mm';

% optional flags
flag_backup = 1; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

pulseq_scanner = 'Siemens_Prisma_3T_Technion';
pns_orientation = 'coronal';

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nxy      = 256;         % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 256  *1e-3;  % [m] FOV geometry
FOV.dz       = 8   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% spiral sequence parameters
SPI_params_mapping();

[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);

SPI.Trec = 5;

%% T2 params
MLEV.n_mlev     = [1 1 1 1 2 2 2 2 4 4 4 4];           % number of MLEV4 preps
MLEV.fSL        = 500;             % [Hz] eff spin-lock field strength
MLEV.t_inter    = 10 *1e-3;         % [s]  inter pulse delay for T2 preparation
MLEV.exc_mode   = 'adiabatic_BIR4'; % 'adiabatic_BIR4' or 'adiabatic_AHP'
MLEV.bir4_tau   = 10 *1e-3;  % [s]  bir4 pulse duration
MLEV.bir4_f1    = 640;       % [Hz] maximum rf peak amplitude
MLEV.bir4_beta  = 10;        % [ ]  am waveform parameter
MLEV.bir4_kappa = atan(10);  % [ ]  fm waveform parameter
MLEV.bir4_dw0   = 30000;     % [rad/s] fm waveform scaling
MLEV = MLEV_init(MLEV, FOV, system);

%% fat saturation and reset
FAT.mode = 'off';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence
for loop_MLEV = 1 : numel(MLEV.n_mlev)   
for loop_NR = 1-SPI.Ndummy : SPI.NR

    % saturtion or crusher
    SAT_add();   
    
    % recovery time
    seq.addBlock(mr.makeDelay(SPI.Trec));

    % fat saturation
    FAT_add();

    % mlev preparation
    MLEV_add();

    % spiral imaging
    SPI_add();
    
end
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();