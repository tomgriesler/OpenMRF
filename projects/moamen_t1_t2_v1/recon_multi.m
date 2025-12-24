%%
clear; clc

%%
study_path = 'E:\University of Michigan Dropbox\Tom Griesler\rawdata_new\technion\2025-12-09\';
study_names = {'meas_MID00493_FID09198_251209_1113_tomgr_cardiac_mrf_3mm.dat'};

%%
for ii = 1:length(study_names)
    study_name_mrf = study_names{ii};
    disp(study_name_mrf)
    clearvars -except study_path study_names study_name_mrf ii
    reco_mrf;
end