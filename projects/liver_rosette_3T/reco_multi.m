%%
clear; clc;

%%
study_path = 'E:\University of Michigan Dropbox\Tom Griesler\rawdata_new\technion\2026-01-28\';
study_names_mrf = dir(fullfile(study_path, 'meas*tomgr_spi_mrf*.dat'));
study_names_mrf = {study_names_mrf.name};

%%
for ii=1:length(study_names_mrf)
    
    clearvars -except ii study_path study_names_mrf

    study_name_mrf = study_names_mrf{ii};

    reco_mrf;
end
