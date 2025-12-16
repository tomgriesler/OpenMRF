%% init pulseq
% basis: SPI readout
% use for: T1 mapping
clear
seq_name = 't1_mapping_5mm';

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

%% inversion
INV.mode = 'on';
% INV.inv_rec_time = [25 45 81 146 263 475 855 1540 2775 5000] *1e-3;
INV.inv_rec_time = [50 100 150 200 250 300 400 500 600 700 800 900 1000 1200 1500 2000 3000 4000 5000] *1e-3;
INV = INV_init(INV, FOV, system);

%% fat saturation and reset
FAT.mode = 'off';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence
for loop_INV = 1 : numel(INV.inv_rec_time)
for loop_NR  = 1-SPI.Ndummy : SPI.NR

    % saturtion or crusher
    SAT_add();   
    
    % recovery time
    seq.addBlock(mr.makeDelay(SPI.Trec));

    % inversion pulse
    INV_add();

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