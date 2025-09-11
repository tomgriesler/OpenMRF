function WASABI = WASABI_init(WASABI, FOV, system)

%% calculate WASABI pulse object and S0 reference delay
WASABI.rf = mr.makeBlockPulse( pi/10, system, ...
                               'Duration', WASABI.tau, ...
                               'freqOffset', 0, ...
                               'PhaseOffset', WASABI.phase, ...
                               'use', 'preparation' );
WASABI.rf.signal = WASABI.rf.signal*0 + WASABI.f1;
WASABI.ref_delay = mr.calcDuration(WASABI.rf);

%% crusher gradients

% calculate crusher objects
if ~isfield(WASABI, 'crush_nTwists_x')
    WASABI.crush_nTwists_x = 4;   % [] number of 2pi twists in x direction
end
if ~isfield(WASABI, 'crush_nTwists_y')
    WASABI.crush_nTwists_y = 4;   % [] number of 2pi twists in y direction
end
if ~isfield(WASABI, 'crush_nTwists_z')
    WASABI.crush_nTwists_z = 11.3;   % [] number of 2pi twists in z direction
end
[WASABI.gx_crush, WASABI.gy_crush, WASABI.gz_crush] = CRUSH_x_y_z(WASABI.crush_nTwists_x, WASABI.crush_nTwists_y, WASABI.crush_nTwists_z, FOV.dx, FOV.dy, FOV.dz, 1/sqrt(3), 1/sqrt(3), system);
WASABI.tcrush = mr.calcDuration(WASABI.gx_crush, WASABI.gy_crush, WASABI.gz_crush);

end

