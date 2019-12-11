% clear; clc, close('all'); warning off

%%
% bs_path = '/opt/matlab_toolboxes/brainstorm3';
% addpath(bs_path);
% brainstorm
% addpath(genpath('./functions'));

close all
set(0,'DefaultFigureWindowStyle','normal')
addpath('/usr/local/MATLAB_Tools/brainstorm3')
brainstorm


anatfile = fullfile(cfg_main.indir,cfg_main.subj,'brainstorm_db/anat', cfg_main.subj);

%%
% datafile = '/data/MEG/Vahab/test_data/raghavan_manoj2/brainstorm_db/anat/raghavan_manoj2';
% datafile = '/data/MEG/Vahab/Scripts/Vahab/Test_MEG2/bednar_peggy/brainstorm_db/anat/bednar_peggy';
% datafile = '/data/MEG/Clinical/MEG/dougherty_danielle/brainstorm_db/anat/dougherty_danielle';
% datafile = '/data/MEG/Vahab/test_data/busby_michael/brainstorm_db/anat/busby_michael';
% datafile = '/data/MEG/Vahab/test_data/krier_henry/brainstorm_db/anat/krier_henry';
% datafile = '/data/MEG/Vahab/test_data/anderson_mark/brainstorm_db/anat/anderson_mark';


disp(anatfile)
d = rdir(fullfile(anatfile,'subjectimage*.mat'));
% d = rdir(fullfile(datafile,'*ras.mat'));
if ~isempty(d)
    sMRI = d.name;
    cd(anatfile)
    cd ..
    OutputFile = fullfile(pwd,'T1.nii');
    out_mri_nii(sMRI, OutputFile);
end

%%
brainstorm stop
rmpath(genpath('/usr/local/MATLAB_Tools/brainstorm3'));
