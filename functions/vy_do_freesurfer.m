function vy_do_freesurfer(cfg_main)

individual_mri = ft_read_mri(cfg_main.mripfile);
ft_sourceplot([], individual_mri);
individual_mri = ft_convert_units(individual_mri, 'mm');

%%
fid = cfg_main.fid;

%%
%         fprintf('Please identify the LPA, RPA, nasion, and a point on the positive Z-axis\n');
%         cfg = [];
%         cfg.method = 'interactive';
%         cfg.coordsys = 'neuromag';
%         individual_mri = ft_volumerealign(cfg, individual_mri);

%% Normalized coordinate system (NCS)
cfg          = [];
% cfg.method   = 'interactive';
cfg.method = 'fiducial'; % the following voxel coords were determined interactive
cfg.coordsys = 'spm';
cfg.fiducial.ac   = fid.NCS.AC;
cfg.fiducial.pc   = fid.NCS.PC;
cfg.fiducial.xzpoint  = fid.NCS.IH;
cfg.fiducial.right   = fid.SCS.RPA;
mri_mni     = ft_volumerealign(cfg, individual_mri);

%%
% reslice & save the transformation matrix to the anatomy_dir
cfg                 = [];
cfg.resolution      = 1;
cfg.dim             = [256 256 256];
mri_resliced        = ft_volumereslice(cfg, mri_mni);

cd(cfg_main.outputmridir)

% filename for saving
mgz_filename = fullfile(cfg_main.outputmridir, [cfg_main.subj, '_mni_resliced' '.mgz']);

% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = mgz_filename;
cfg.filetype        = 'mgz';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_resliced)

% Save the transformation matrix
transform_vox2mni   = mri_resliced.transform;
filename_vox2mni    = fullfile(cfg_main.outputmridir, [cfg_main.subj, '_transform_vox2mni']);
save(filename_vox2mni, 'transform_vox2mni');

%% Subject Coordinate System (neuromag)
cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'neuromag';
cfg.fiducial.nas   = fid.SCS.NAS;
cfg.fiducial.lpa   = fid.SCS.LPA;
cfg.fiducial.rpa   = fid.SCS.RPA;
cfg.spmversion     = 'spm12';
mri_neuromag = ft_volumerealign(cfg, mri_resliced);

%%
headshape = ft_read_headshape(cfg_main.hsfile);
headshape = ft_convert_units(headshape, 'mm');

cfg = [];
cfg.method = 'headshape';
cfg.headshape.interactive = 'no';
cfg.headshape.icp = 'yes';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'neuromag';
cfg.spmversion     = 'spm12';
mri_neuromag = ft_volumerealign(cfg, mri_neuromag);

%% == to check everything went right with co-reg!
ft_determine_coordsys(mri_neuromag, 'interactive', 'no')
ft_plot_headshape(headshape);
view([-90, 0]),

%%
% Save the transformation matrix
transform_vox2neuromag   = mri_neuromag.transform;
filename_vox2neuromag    = fullfile(cfg_main.outputmridir, [cfg_main.subj, '_transform_vox2neuromag']);
save(filename_vox2neuromag, 'transform_vox2neuromag');

%%
% FSL variables
threshold       = 0.5;
T               = inv(mri_resliced.transform);
center          = round(T(1:3,4))';
subjectname     = cfg_main.subj;

% name for the temporary nifti file
t   = fullfile(cfg_main.outputmridir, [cfg_main.subj, '_nifti_tmp']);

%%
anatomy_dir = cfg_main.outputmridir;
subject_code = cfg_main.subj;

% Convert to nifti temporarily and save;
cfg = [];
cfg.filename = t;
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_resliced);

% Create the FSL command-string
str = ['/opt/fsl/bin/bet ',t,'.nii ',t];
str = [str,'-R -f ',num2str(threshold),' -c ', num2str(center),' -g 0 -m -v'];

% Call the FSL command-string
system(str);

% Read the FSL-based segmentation
seg  = ft_read_mri([t,'-R.nii.gz']);
delete([t,'.nii']);
delete([t,'-R.nii.gz']);
delete([t,'-R_mask.nii.gz']);

mri_skullstrip               = fullfile(anatomy_dir, [subject_code, '_skullstrip']);

% Save the FSL-based segmentation in .mgz
cfg = [];
cfg.filename = mri_skullstrip;
cfg.filetype = 'mgz';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, seg);

% Check the plot already now
skullstrip = ft_read_mri([mri_skullstrip '.mgz']);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, skullstrip);


%% Freesurfer scripts (creates subject-specific subdirectory in the directory where previous files are stored)

ft_path = ('/data/MEG/Vahab/Github/fieldtrip');
if ~ft_hastoolbox('qsub',1)
    addpath(fullfile(ft_path,'qsub'));
end

% subjects = {'s18' 's15' 's27' 's28'};
%
% for i = 1:numel(subjects)
%
%   subject = subjects{i};
%
%   qsubfeval('qsub_streams_anatomy_freesurfer', subject,...
%             'memreq', 1024^3 * 6,...
%             'timreq', 720*60,...
%             'batchid', 'streams_freesurferI');
% end

subject = cfg_main.subj;
% Fressurfer script1
shell_script      = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/vy_anatomy_freesurfer.sh';
mri_dir           = anatomy_dir;
subject_dir       = subject;

% create the string pointing to streams_anatomy_freesurfer.sh
command = [shell_script, ' ', mri_dir, ' ', subject_dir];

% call the script
system(command);

% running under command
% recon-all -s $subject -sd $SUBJECTS_DIR -make all
% recon-all -s bednar_peggy -sd /data/MEG/Clinical/ft_process/19/bednar_peggy/anat -make all

%% Check-up and white matter segmentation cleaning if needed
anatomy_dir     = '/data/MEG/Clinical/ft_process/19/bednar_peggy/anat/BAK';
inp_dir         = fullfile(anatomy_dir, subject_dir);

t1              = fullfile(inp_dir, 'mri', 'T1.mgz');
normalization2  = fullfile(inp_dir, 'mri', 'brain.mgz');
white_matter    = fullfile(inp_dir, 'mri', 'wm.mgz');
white_matter_old = fullfile(inp_dir, 'mri', 'wm_old.mgz');

% Show T1
mri = ft_read_mri(t1);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subject.name ' ' 'T1'], 'numbertitle', 'off');

% Show skullstripped image
mri = ft_read_mri(normalization2);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subject.name ' ' 'skull-stripped'], 'numbertitle', 'off');

% Show white matter image
mri = ft_read_mri(white_matter);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subject.name ' ' 'white matter'], 'numbertitle', 'off');

if exist(white_matter_old)
    
    mri = ft_read_mri(white_matter_old);
    cfg = [];
    cfg.interactive = 'yes';
    ft_sourceplot(cfg, mri);
    set(gcf, 'name', [subject.name ' ' 'white matter old'], 'numbertitle', 'off');
    
end

%%
cfg = [];
cfg.mri_dir = anatomy_dir;
vy_anatomy_wmclean(cfg, subject)

% streams_anatomy_volumetricQC(subject)
%
% streams_anatomy_wmclean(subject)

%%
subject = cfg_main.subj;
% Fressurfer script1
shell_script      = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/vy_anatomy_freesurfer2.sh';
mri_dir           = anatomy_dir;
subject_dir       = subject;

% create the string pointing to streams_anatomy_freesurfer.sh
command = [shell_script, ' ', mri_dir, ' ', subject_dir];

% call the script
system(command);

%% Post-processing Freesurfer script: workbench HCP tool
ft_path = ('/data/MEG/Vahab/Github/fieldtrip');
if ~ft_hastoolbox('qsub',1)
    addpath(fullfile(ft_path,'qsub'));
end
addpath(genpath(fullfile('/opt/workbench/')));

% Strings for the command
shell_script      = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/vy_anatomy_postfreesurferscript.sh';
mri_dir           = anatomy_dir;
subject_dir       = subject;

% streams_anatomy_freesurfer2.sh
command = [shell_script, ' ', mri_dir, ' ', subject_dir];

system(command);

%%  Sourcemodel
cfg = [];
cfg.anatomy_dir = mri_dir;
vy_anatomy_sourcemodel2d(cfg,subject)

%% Headmodel
cfg = [];
cfg.anatomy_savedir = mri_dir;
vy_anatomy_headmodel(cfg,subject);

%%  Coregistration check
cfg = [];
cfg.anatomy_dir = mri_dir;
vy_anatomy_coregistration_qc(cfg, subject);

%% Leadfield parcellation
cfg = [];
cfg.anatomy_dir = mri_dir;
cfg.megdata = cfg_main.megdata;
vy_leadfield(cfg,subject)

%%
load(fullfile(anatomy_dir,[subject,'_headmodel.mat']));
load(fullfile(anatomy_dir,[subject,'_leadfield.mat']));
load(fullfile(anatomy_dir,[subject,'_sourcemodel.mat']));

figure;
ft_plot_vol(headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
hold on;
ft_plot_headshape(headshape);
ft_plot_mesh(leadfield.pos(leadfield.inside, :));
view ([0 90])
