
%% 
% ft_path = ('/opt/matlab_toolboxes/fieldtrip');
% addpath(ft_path);
% ft_defaults

%%
cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, ~] = vy_init(cfg_init);

%% Read the DICOM files
% dicomfile = '/data/MEG/Clinical/MRI/xxx/DICOM/EXP00000/EXP0000';
[f, p] = uigetfile('*');
dicomfile = fullfile(p, f);
mri = ft_read_mri(dicomfile);

%%
ft_sourceplot([], mri);
mri = ft_convert_units(mri, 'mm');

%% Filename for saving
cd(p)
cd ..
cd ..
nii_filename = 'mri.nii';

%% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = nii_filename;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri);

%%
addpath(genpath(allpath.spm_path))
spm_get_defaults
setorigin_center(nii_filename);

%%
mri1 = ft_read_mri(nii_filename);
ft_sourceplot([], mri1);

%% Reslice & save the transformation matrix to the anatomy_dir
cfg                 = [];
cfg.resolution      = 1;
cfg.dim             = [256 256 256];
mri_resliced        = ft_volumereslice(cfg, mri1);
ft_sourceplot([], mri_resliced);

%%
cfg          = [];
cfg.method   = 'interactive';
cfg.coordsys = 'spm';
mri_mni     = ft_volumerealign(cfg, mri_resliced);

%% Filename for saving
mgz_filename = 'mni_resliced.mgz';

%% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = mgz_filename;
cfg.filetype        = 'mgz';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_resliced);

%% Save the transformation matrix
transform_vox2mni   = mri_resliced.transform;
filename_vox2mni    = 'transform_vox2mni';
save(filename_vox2mni, 'transform_vox2mni');

%%
