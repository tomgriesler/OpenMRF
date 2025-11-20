%% init pulseq
% basis: SPI readout
% use for: T1 mapping
clear
seq_name = 't1_mapping_tse';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

pulseq_scanner = 'Siemens_Prisma_3T_Technion';

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nx       = 256;
FOV.Ny       = 256;
FOV.fov_x    = 256 *1e-3;
FOV.fov_y    = 256 *1e-3;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% TSE sequence parameters
TSE.Nrep      = 1;
TSE.n_echo    = 1;
TSE.TE        = 11    *1e-3;
TSE.TR        = 0  *1e-3;
TSE.Ndummy    = 4;
TSE.exc_time  = 3.0 *1e-3;
TSE.rfc_time  = 3.0 *1e-3;
TSE.exc_tbw   = 4;
TSE.rfc_tbw   = 4;
TSE.t_acq     = 6.4 *1e-3 + 2*system.adcDeadTime;
TSE.os_mode   = 1;  % read oversampling: 0 off, 1 on
TSE.mode_exc  = 'sinc'; % 'sigpy_SLR' or 'sinc'
TSE.mode_rfc  = 'sinc'; % 'sigpy_SLR' or 'sinc'
TSE.enc_mode  = 'centric';

[TSE, ktraj_adc, ktraj_full] = TSE_init(TSE, FOV, system);

TSE.Trec = 2;

%% inversion
INV.mode = 'on';
INV.inv_rec_time = [25] *1e-3; %  45 81 146 263 475 855 1540 2775 5000] *1e-3; % 
INV = INV_init(INV, FOV, system);

%% fat saturation and reset
FAT.mode = 'off';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence
for loop_INV = 1 : numel(INV.inv_rec_time)
for loop_TR  = 1-TSE.Ndummy : TSE.nex

    % saturtion or crusher
    SAT_add();   
    
    % recovery time
    seq.addBlock(mr.makeDelay(TSE.Trec));

    % inversion pulse
    INV_add();

    % fat saturation
    FAT_add();

    % spiral imaging
    TSE_add();  

end    
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();