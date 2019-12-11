%% Check-up and white matter segmentation cleaning if needed
% inp_dir         = '/data/MEG/Clinical/MRI/bednar_peggy/bednar_peggy_04_04_2019_FSRecon';
inp_dir         = '/MEG_data/Vahab/Processed_data/ft_process/corcoran_margaret/anat/corcoran_margaret';
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
% ft_path = ('/data/MEG/Vahab/Github/fieldtrip');
if ~ft_hastoolbox('qsub',1)
    addpath(fullfile(ft_path,'qsub'));
end
addpath(genpath(fullfile('/opt/workbench/')));

% Strings for the command
% shell_script      = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/vy_anatomy_postfreesurferscript.sh';
% mri_dir           = '/data/MEG/Clinical/ft_process/19/bednar_peggy/anat/test';

shell_script      = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/functions/vy_anatomy_postfreesurferscript_megneto.sh';
mri_dir           = '/MEG_data/Vahab/Processed_data/ft_process/corcoran_margaret/anat';
subject_dir       = subject.name;

% streams_anatomy_freesurfer2.sh
command = [shell_script, ' ', mri_dir, ' ', subject_dir];

system(command);


%%
mgz_filename = fullfile(inp_dir,'mri/T1.mgz');
mri_native = ft_read_mri(mgz_filename);

nii_filename1 = fullfile(mri_dir, [subj, '_native' '.nii']);
cfg                 = [];
cfg.filename        = nii_filename1;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_native);

%%
nii_filename = fullfile(mridir,'T1.nii');
mri_bs = ft_read_mri(nii_filename);
ft_sourceplot([], mri_bs);

nii_filename2 = fullfile(mri_dir, [subj, '_bs' '.nii']);
cfg                 = [];
cfg.filename        = nii_filename2;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_bs);

%% Co-reg, estimate and reslice
addpath(genpath(allpath.spm_path))
spm_get_defaults

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[nii_filename2,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[nii_filename1,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);

%%
rmpath(genpath(allpath.spm_path))

%%
nii_filename3 = fullfile(mri_dir, ['r',subj, '_native' '.nii']);
rmri_native = ft_read_mri(nii_filename3);
ft_sourceplot([], rmri_native);
% cfg                 = [];
% cfg.resolution      = 1;
% cfg.dim             = [256 256 256];
% mri_resliced        = ft_volumereslice(cfg, rmri_native);
% ft_sourceplot([], mri_resliced);
mri_resliced    =  rmri_native;

%%
cfg          = [];
% cfg.method   = 'interactive';
cfg.method = 'fiducial'; % the following voxel coords were determined interactive
cfg.coordsys = 'spm';
cfg.fiducial.ac   = fid.NCS.AC;
cfg.fiducial.pc   = fid.NCS.PC;
cfg.fiducial.xzpoint  = fid.NCS.IH;
cfg.fiducial.right   = fid.SCS.RPA;
mri_native_tr     = ft_volumerealign(cfg, mri_resliced);
% ft_sourceplot([], mri_native_tr);

%%
cfg = [];
cfg.method = 'fiducial';
% cfg.method   = 'interactive';
cfg.coordsys = 'neuromag';
cfg.fiducial.nas   = fid.SCS.NAS;
cfg.fiducial.lpa   = fid.SCS.LPA;
cfg.fiducial.rpa   = fid.SCS.RPA;
cfg.spmversion     = 'spm12';
mri_native_tr = ft_volumerealign(cfg, mri_native_tr);
% ft_sourceplot([], individual_mri_mni_neuromag);

%%
cfg = [];
cfg.method = 'headshape';
cfg.headshape.interactive = 'no';
cfg.headshape.icp = 'yes';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'spm';
cfg.spmversion     = 'spm12';
mri_realigned = ft_volumerealign(cfg, mri_native_tr);
% ft_sourceplot([], mri_realigned);

%% To check everything went right with co-reg!
ft_determine_coordsys(mri_realigned, 'interactive', 'no')
ft_plot_headshape(headshape);
view([-90, 0]),

%%
transform_vox   = mri_resliced.transform;
filename_vox    = fullfile(mri_dir, [subj, '_transform_vox']);
save(filename_vox, 'transform_vox');

transform_vox2spm   = mri_realigned.transform;
filename_vox2spm    = fullfile(mri_dir, [subj, '_transform_vox2spm']);
save(filename_vox2spm, 'transform_vox2spm');

%%  Sourcemodel
cfg = [];
cfg.anatomy_dir = mri_dir;
sourcemodel = vy_anatomy_sourcemodel2d_test(cfg, subj);

%% Headmodel
cfg = [];
cfg.anatomy_savedir = mri_dir;
headmodel = vy_anatomy_headmodel(cfg, subj);

%%  Coregistration check
cfg = [];
cfg.anatomy_dir = mri_dir;
vy_anatomy_coregistration_qc(cfg, subj);

%% Leadfield parcellation
cfg = [];
% cfg.atlasdir = '/data/MEG/Vahab/Github/MCW-MEGlab/tools/dyncon_erfosc';
cfg.anatomy_dir = mri_dir;
cfg.megdata = t_data.all;
vy_leadfield(cfg,subj)

%%

load(fullfile(mri_dir,[subj,'_headmodel.mat']));
load(fullfile(mri_dir,[subj,'_leadfield.mat']));
load(fullfile(mri_dir,[subj,'_sourcemodel.mat']));

T1 = mri_resliced.transform;
T2 = mri_realigned.transform;

sourcemodel1 = ft_transform_geometry(T1/mri_native.transform, sourcemodel);
sourcemodel2 = ft_transform_geometry((T2/T1), sourcemodel1);

leadfield1 = ft_transform_geometry((T2/T1), leadfield);
headmodel1 = ft_transform_geometry((T2/T1), headmodel);

figure; hold on;
ft_plot_vol(headmodel1, 'facecolor', 'none'); alpha 0.5;
ft_plot_mesh(sourcemodel2, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
ft_plot_mesh(leadfield1.pos(leadfield1.inside, :));
view ([0 90])
hold on;
ft_plot_headshape(headshape);

%%
cfg = [];
cfg.headmodel = headmodel1;
cfg.sourcemodel = sourcemodel2;
cfg.grid   = leadfield1;
cfg.mtag = 'dics_fs';
cfg.sens = ep_data.all.grad;
cfg.subj = subj;
cfg.outputdir = outd.sub;
vy_source_dics_fs(cfg, ep_data);

%%
% load(fullfile(mri_dir,[subj,'_headmodel.mat']));
% load(fullfile(mri_dir,[subj,'_leadfield.mat']));
% load(fullfile(mri_dir,[subj,'_sourcemodel.mat']));
% 
% headshape1 = ft_transform_geometry(transform_vox2neuromag/transform_vox,headshape);
% 
% figure;
% ft_plot_vol(headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% hold on;
% ft_plot_headshape(headshape1);
% ft_plot_mesh(leadfield.pos(leadfield.inside, :));
% view ([0 90])
% ft_plot_headshape(sourcemodel2);


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
% cfg = [];
% cfg.anatomy_dir = mri_dir;
% vy_anatomy_sourcemodel2d(cfg,subject.name);

%%


