%% add magnetization inversion pulse
if ~exist('loop_INV', 'var')
    loop_INV = 1;
end
seq.addBlock(INV.rf);
seq.addBlock(INV.gx_crush, INV.gy_crush, INV.gz_crush);
seq.addBlock(INV.inv_rec_delay(loop_INV));