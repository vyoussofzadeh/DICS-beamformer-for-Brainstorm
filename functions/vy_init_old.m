function allpath = vy_init(cfg)

cd_org = cd;
addpath(genpath(cd_org));

%- fieldtrip
% ft_path = fullfile(cfg.path_tools,'/ft packages/fieldtrip-master');
ft_path = fullfile(cfg.path_tools,'/ft_packages/fieldtrip_20190419');

addpath(ft_path);
ft_defaults

ft_old = fullfile(cfg.path_tools,'ft_packages/ft_old');
%- colormap
addpath(fullfile(ft_path, 'external/brewermap'));
%-fastICA
addpath(ft_path,'/external/fastica');

%- HCP
hcp_path = fullfile(cfg.path_tools,'/megconnectome-3.0');
addpath(genpath(hcp_path));

%- atlas
atlas_path = fullfile(ft_path,'template','atlas');
atlas = ft_read_atlas(fullfile(atlas_path,'aal/ROI_MNI_V4.nii'));

%-Grid template
load temp_grid_8mm % from, vy_warping()
% template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %
% template_mri = ft_read_mri(fullfile(ft_path,'template/anatomy','single_subj_T1.nii')); %

%-CONN
connpath = fullfile(cfg.path_tools,'/tools/Conn/conn');
% addpath(connpath);

%-SPM
% spm_path = '/opt/matlab_toolboxes/spm12';
spm_path = fullfile(cfg.path_tools,'SPM/spm12');
% addpath(spm_path);

%-SPM-beamformer
spmbf_path = fullfile(cfg.path_tools,'Beamforming');

%%
allpath.cd_org = cd_org;
allpath.ft_path = ft_path;
allpath.ft_old = ft_old;
allpath.hcp_path = hcp_path;
allpath.atlas_path = atlas_path;
% allpath.template_mri = template_mri;
allpath.connpath = connpath;
allpath.spm_path = spm_path;
allpath.spmbf_path = spmbf_path;



