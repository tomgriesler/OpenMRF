%% init
% basis: spiral readout
% use for: B0 and B1 mapping via WASABI method -> doi.org/10.1002/mrm.26133
clear
seq_name = 'spi_wasabi_5mm';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 1; % 0: off,  1: simulate PNS stimulation
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
FOV.dz       = 5  *1e-3;    % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% spiral sequence parameters
SPI_params_wasabi();

SPI.lim_reph_slew = 0.5;        % [ ] reduce stimulation during gradient spoiling

[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);
SPI.Trec = 0.5;

%% fat saturation and reset
FAT.mode = 'off';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% WASABI preparation
WASABI.f0       = 123216135;     % [Hz]  larmor frequency
WASABI.ppm      = -2 : 0.1 : 2;  % [ppm] wasabi pulse offresonances
WASABI.f1       = 150;           % [Hz]  wasabi pulse amplitude
WASABI.tau      = 5 *1e-3;       % [s]   pulse duration (5ms)
WASABI.phase    = 0;             % [rad] pulse phase
WASABI.B0       = WASABI.f0 / system.gamma; % [T]
WASABI.B1       = WASABI.f1 / system.gamma; % [T]
WASABI.n_ppm    = numel(WASABI.ppm);        % [ ]
WASABI.f_off    = WASABI.ppm *1e-6 *WASABI.B0 *system.gamma; % [Hz] wasabi pulse offresonace frequencies
WASABI.tau      = round(WASABI.tau/system.rfRasterTime) * system.rfRasterTime; % prevent timing errors
WASABI          = WASABI_init(WASABI, FOV, system);

%% add sequence loops

ndummy = 1;
for loop_ppm = 1 : WASABI.n_ppm + 1
for loop_NR = 1-ndummy:SPI.NR

    ndummy = 0;

    % saturation
    SAT_add();

    % recovery time
    seq.addBlock(mr.makeDelay(SPI.Trec));

    % fat saturation
    FAT_add();

    % WASABI preparation
    WASABI_add();

    % spiral readouts
    SPI_add();

end    
end

%% plot sequence diagram
seq.plot('TimeRange',[2 5]*SPI.Trec)

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();