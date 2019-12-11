%% Check-up and white matter segmentation cleaning if needed
inp_dir         = '/data/MEG/Clinical/MRI/bednar_peggy/bednar_peggy_04_04_2019_FSRecon';
subject.name = subj;

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
% anatomy_dir = '/data/MEG/Clinical/ft_process/19/bednar_peggy/anat/test/bednar_peggy_04_04_2019_FSRecon';
% 
% cfg = [];
% cfg.mri_dir = anatomy_dir;
% vy_anatomy_wmclean_test(cfg, subject)

%% Post-processing Freesurfer script: workbench HCP tool
ft_path = ('/data/MEG/Vahab/Github/fieldtrip');
if ~ft_hastoolbox('qsub',1)
    addpath(fullfile(ft_path,'qsub'));
end
addpath(genpath(fullfile('/opt/workbench/')));

% Strings for the command
shell_script      = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/vy_anatomy_postfreesurferscript.sh';
mri_dir           = '/data/MEG/Clinical/ft_process/19/bednar_peggy/anat/test';
subject_dir       = subject.name;

% streams_anatomy_freesurfer2.sh
command = [shell_script, ' ', mri_dir, ' ', subject_dir];

system(command);

%% Native space
% filename for saving
mgz_filename = '/data/MEG/Clinical/ft_process/19/bednar_peggy/anat/test/bednar_peggy/mri/T1.mgz';
individual_mri = ft_read_mri(mgz_filename);

% Save the transformation matrix
transform_vox   = individual_mri.transform;
filename_vox    = fullfile(outputmridir, [subject_dir, '_transform_vox']);
save(filename_vox, 'transform_vox');

% %% mni, neuromag space
% nii_filename = '/data/MEG/Clinical/MRI/bednar_peggy/T1.nii';
% individual_mri_bs = ft_read_mri(nii_filename);

% filename for saving
mgz_filename = fullfile(mri_dir, [subj, '_mni_resliced' '.mgz']);

% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = mgz_filename;
cfg.filetype        = 'mgz';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, individual_mri)

%%
individual_mri = ft_read_mri(cfg_main.mripfile);
%     ft_sourceplot([], individual_mri);
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
cfg.coordsys = 'neuromag';
cfg.fiducial.ac   = fid.NCS.AC;
cfg.fiducial.pc   = fid.NCS.PC;
cfg.fiducial.xzpoint  = fid.NCS.IH;
cfg.fiducial.right   = fid.SCS.RPA;
mri_realigned     = ft_volumerealign(cfg, individual_mri);

%% Subject Coordinate System (SCS / CTF)
cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'neuromag';
cfg.fiducial.nas   = fid.SCS.NAS;
cfg.fiducial.lpa   = fid.SCS.LPA;
cfg.fiducial.rpa   = fid.SCS.RPA;
cfg.spmversion     = 'spm12';
mri_realigned = ft_volumerealign(cfg, mri_realigned);

%%
%         ft_sourceplot([], mri_realigned); hold on
%         plot3(fid.SCS.NAS(1,1), fid.SCS.NAS(1,2), fid.SCS.NAS(1,3), 'm*');
%         plot3(fid.SCS.LPA(1,1), fid.SCS.LPA(1,2), fid.SCS.LPA(1,3), 'm*');
%         plot3(fid.SCS.RPA(1,1), fid.SCS.RPA(1,2), fid.SCS.RPA(1,3), 'm*');

%%
headshape = ft_read_headshape(datafile);
headshape = ft_convert_units(headshape, 'mm');

%%
cfg = [];
cfg.method = 'headshape';
cfg.headshape.interactive = 'no';
cfg.headshape.icp = 'yes';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'neuromag';
cfg.spmversion     = 'spm12';
mri_realigned = ft_volumerealign(cfg, mri_realigned);

%%
cfg = [];
cfg.method = 'headshape';
cfg.headshape.interactive = 'no';
cfg.headshape.icp = 'yes';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'neuromag';
cfg.spmversion     = 'spm12';
mri_realigned = ft_volumerealign(cfg, mri_realigned);

%% == to check everything went right with co-reg!
ft_determine_coordsys(mri_realigned, 'interactive', 'no')
ft_plot_headshape(headshape);
view([-90, 0]),

%%
% Save the transformation matrix
transform_vox2neuromag   = individual_mri.transform;
filename_vox2neuromag    = fullfile(cfg_main.outputmridir, [cfg_main.subj, '_transform_vox2neuromag']);
save(filename_vox2neuromag, 'transform_vox2neuromag');

%%  Sourcemodel
cfg = [];
cfg.anatomy_dir = mri_dir;
vy_anatomy_sourcemodel2d_test(cfg,subj)

%% Headmodel
cfg = [];
cfg.anatomy_savedir = mri_dir;
vy_anatomy_headmodel(cfg, subj);

%%  Coregistration check
cfg = [];
cfg.anatomy_dir = mri_dir;
vy_anatomy_coregistration_qc(cfg, subj);

%%
%% Leadfield parcellation
cfg = [];
cfg.anatomy_dir = mri_dir;
cfg.megdata = t_data.all;
vy_leadfield(cfg,subj)

%%
load(fullfile(mri_dir,[subj,'_headmodel.mat']));
load(fullfile(mri_dir,[subj,'_leadfield.mat']));
load(fullfile(mri_dir,[subj,'_sourcemodel.mat']));

figure;
ft_plot_vol(headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
hold on;
ft_plot_headshape(headshape);
ft_plot_mesh(leadfield.pos(leadfield.inside, :));
view ([0 90])

%%
% %% Normalized coordinate system (NCS)
% cfg          = [];
% % cfg.method   = 'interactive';
% cfg.method = 'fiducial'; % the following voxel coords were determined interactive
% cfg.coordsys = 'neuromag';
% cfg.fiducial.ac   = fid.NCS.AC;
% cfg.fiducial.pc   = fid.NCS.PC;
% cfg.fiducial.xzpoint  = fid.NCS.IH;
% cfg.fiducial.right   = fid.SCS.RPA;
% mri_realigned     = ft_volumerealign(cfg, individual_mri);
%
% %% Subject Coordinate System (SCS / CTF)
% cfg = [];
% cfg.method = 'fiducial';
% cfg.coordsys = 'neuromag';
% cfg.fiducial.nas   = fid.SCS.NAS;
% cfg.fiducial.lpa   = fid.SCS.LPA;
% cfg.fiducial.rpa   = fid.SCS.RPA;
% cfg.spmversion     = 'spm12';
% mri_realigned = ft_volumerealign(cfg, mri_realigned);
%
%
% %%
% headshape = ft_read_headshape(datafile);
% headshape = ft_convert_units(headshape, 'mm');
%
% cfg = [];
% cfg.method = 'headshape';
% cfg.headshape.interactive = 'no';
% cfg.headshape.icp = 'yes';
% cfg.headshape.headshape = headshape;
% cfg.coordsys = 'neuromag';
% cfg.spmversion     = 'spm12';
% mri_neuromag = ft_volumerealign(cfg, mri_realigned);
%
% %% == to check everything went right with co-reg!
% ft_determine_coordsys(mri_neuromag, 'interactive', 'no')
% ft_plot_headshape(headshape);
% view([-90, 0]),

%%  Sourcemodel
cfg = [];
cfg.anatomy_dir = mri_dir;
vy_anatomy_sourcemodel2d(cfg,subject.name);

%%


