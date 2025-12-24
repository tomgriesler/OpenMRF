%% init pulseq
% basis: SPI readout
% use for: T2 mapping
clear
seq_name = 't2_mapping_5mm';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

pulseq_scanner = 'Siemens_Prisma_3T_Technion';
pns_orientation = 'coronal';

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nxy      = 96;         % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 96  *1e-3;  % [m] FOV geometry
FOV.dz       = 5   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% spiral sequence parameters
SPI_params_mapping();

SPI.lim_spoil_slew = 0.6;        % [ ] reduce stimulation during gradient spoiling

[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);

SPI.Trec = 5;

%% T2 params
T2.exc_mode   = 'adiabatic_BIR4'; 
T2.rfc_dur    = 2 *1e-3;
% T2.prep_times = [16	32	48	64	80	96	112	128	144	160] * 1e-3;
% T2.prep_times = [12.6 25.2 37.8 50.4 63 75.6 88.2 100.8 113.4 126 138 ...
%     151.2 163.8 176.4 189 201.6 214.2 226.8 239.4 252 264.6 277.2 289.8 ...
%     302.4 315 327.6 340.2 352.8 365.4 378 390.6 403.2] * 1e-3;
% T2.prep_times = (25.2:25.2:403.2) * 1e-3;
T2.prep_times = (10.9:10.9:349) * 1e-3;
T2            = T2_init(T2, FOV, system);

%% fat saturation and reset
FAT.mode = 'off';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence
for loop_T2 = 1 : numel(T2.prep_times)   
for loop_NR = 1-SPI.Ndummy : SPI.NR

    % saturtion or crusher
    SAT_add();   
    
    % recovery time
    seq.addBlock(mr.makeDelay(SPI.Trec));

    % fat saturation
    FAT_add();

    % mlev preparation
    T2_add();

    % spiral imaging
    SPI_add();
    
end
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();