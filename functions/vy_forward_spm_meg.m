function vy_forward_spm_meg(cfg_main)

val = 1;

%% initial settings
restoredefaultpath
% addpath(genpath(spm_path))
addpath(genpath(cfg_main.allpath.spm_path))
spm_get_defaults
addpath(genpath(cfg_main.allpath.spmbf_path))
addpath(genpath(cfg_main.allpath.hcp_path))
% addpath(genpath(cfg_main.p.ft_path))

%%
% d = ['.\m',subj];
% if exist([d,'.mat'], 'file') == 2
%     D = spm_eeg_load(d);
% else
S = [];
S.dataset = cfg_main.datafile;
S.channels = {'meg'};
D1 = spm_eeg_convert(S);
% spm_eeg_review(D1)
D1.save
% end

%%
D = spm_eeg_ft2spm(cfg_main.eint_data, cfg_main.sub);

D2 = struct(D);
D2.fiducials = D1.fiducials;
D = meeg(D2);
% spm_eeg_review(D)

%% mri segmentation & mesh generation
% sMRI = fullfile(mridir,[subj,'_T1.nii']);
sMRI = cfg_main.mripfile;
% if exist(sMRI, 'file') ~= 2
%     sMRI = 'F:\My Matlab\SPM\spm12_4\spm12\spm12\canonical\single_subj_T1.nii';
% end
Msize = 2; val = 1;
D = spm_eeg_inv_mesh_ui(D, val, sMRI, Msize);
D.save

%%
selection = 1:3;
meegfid = D1.fiducials;
meegfid.fid.pnt   = meegfid.fid.pnt(selection, :);
meegfid.fid.label = meegfid.fid.label(selection);
newmrifid = D.inv{val}.mesh.fid;
newmrifid.fid.pnt   = D.inv{val}.mesh.fid.fid.pnt(1:3,:);
newmrifid.fid.label = {'Nasion';'LPA';'RPA'}; %{'Nasion';'LPA';'RPA'};
useheadshape = 1; val = 1;
D = spm_eeg_inv_datareg_ui(D, val, meegfid, newmrifid, useheadshape);
% D = spm_eeg_inv_datareg_ui(D);

%% forward modelling
D.inv{1}.forward(1).voltype = 'Single Shell';
D = spm_eeg_inv_forward(D);
spm_eeg_inv_checkforward(D,1);
% D = spm_eeg_inv_forward_ui(D);
D.save

savepath = ['ForwardModel_',cfg_main.sub];
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

%% source BF
S = [];
S.D = D;
S.trigger{1} = {'Undefined'};
S.timewindows{1}=[0.3;0.7];
S.timewindows{2}=[-0.4;0];
S.trigger{2} = {'Undefined'};
% S.freqbands{1} = [18,23];
S.freqbands{1} = [2,23];
% S.freqbands{1} = [1,8];
S.gridstep = 10;
S.regpc = 5;
S.preview = 1;
spm_eeg_ft_beamformer_gui(S);

%% Avergaring
% datafile3 = ['.\aefspmeeg_',name];
S = [];
S.D = D;
S.robust = false;
S.prefix = 'm';
D = spm_eeg_average(S);
D.save
% spm_eeg_review(D)

% end


%% source task
% inv_typ = 'GS';% 'GS'; %'COH';   % (MSP)%inv_typ = 'IID';    % Minimum Norm Least-Squares
inv_typ = 'COH'; % Loreta
woi = [0 999];
D.inv{val}.inverse = [];
D.inv{val}.inverse.type   = inv_typ;
D.inv{val}.inverse.lpf    = 0;
D.inv{val}.inverse.hpf    = 48;
D.inv{val}.inverse.woi    = woi;
D.inv{val}.inverse.Han    = 1;
D.inv{val}.inverse.modality = {'meg'};
D = spm_eeg_invert(D);
D.save

% con_win = {[400,700];0};
% woi1 = [0,300];
% con_win = {woi1; 0};
% D.inv{val}.contrast.woi  = con_win{1};
% D.inv{val}.contrast.fboi = con_win{2};
% D.inv{val}.contrast.type = 'evoked';
% D = spm_eeg_inv_results(D);
%
% savepath = ['t',num2str(con_win{1}(1)),'_',num2str(con_win{1}(2))];
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
%
% D.inv{val}.contrast.smooth = 12;
% D.inv{val}.contrast.format = {'image'};
% D = spm_eeg_inv_Mesh2Voxels(D);
% D.save

% con_win = {[400,700];0};
woi1 = [400,700];
con_win = {woi1; 0};
D.inv{val}.contrast.woi  = con_win{1};
D.inv{val}.contrast.fboi = con_win{2};
D.inv{val}.contrast.type = 'evoked';
D = spm_eeg_inv_results(D);

savepath = ['t',num2str(con_win{1}(1)),'_',num2str(con_win{1}(2))];
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

D.inv{val}.contrast.smooth = 12;
D.inv{val}.contrast.format = {'image'};
D = spm_eeg_inv_Mesh2Voxels(D);
D.save
% D = spm_eeg_invert_ui(D);
% spm_eeg_inv_imag_api(D)


%%
