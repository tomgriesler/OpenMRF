% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ----------------- readout: spiral (SPI) -----------------
% -------------------- Noise pre-scans --------------------
% ---------------------------------------------------------

if isfield(SPI, 'Nnoise')
    for loop_noise = 1:SPI.Nnoise
        [seq, TRID] = GE_add_TRID(seq, TRID, 'noise_prescan', flag_GE);
        seq.addBlock(SPI.adc, ceil((mr.calcDuration(SPI.adc) + 1e-3) / system.blockDurationRaster) * system.blockDurationRaster);
    end
end
