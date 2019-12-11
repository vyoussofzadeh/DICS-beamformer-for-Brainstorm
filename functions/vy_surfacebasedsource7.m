clear fid
if ~isempty(d)
    sMRI1 = d.name;
    load(sMRI1);
    fid.SCS = SCS;
    fid.NCS = NCS;
    mripfile = fullfile(mridir,'T1.nii');
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    
    cfg = [];
    cfg.megdata = t_data.app;
    cfg.mripfile = mripfile;
    cfg.hsfile = datafile; % headshape;
    cfg.fid = fid;
    cfg.outputmridir = outputmridir;
    cfg.subj = subj;
    cfg.plotflag = 2;
    outanat = vy_mri_neuromag2(cfg);
    
end


%%
bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
bsanatdir = fullfile(indir,subj,'brainstorm_db/anat');
fsdir = '/MEG_data/MRI_database/epilepsy';

%%
anatomy_dir = '/data/MEG/Vahab/test_data/fs_test'; % Squiggle server
subject.name = subj;
mri_dir      = bsanatdir;
cd(fullfile(mri_dir))

%%
d = rdir(fullfile(fsdir,['*',subj(end-5:end)],'*_FSRecon/mri/T1.mgz'));
length(d)
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    mgz_filename{i} = d(i).name;
end
disp(mgz_filename')


%%
scan = input('choose scan:');
nsub = length(subsel);
disp('============');

%%
mgz = mgz_filename(scan);
mri_native = ft_read_mri(mgz{1});

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

%%
% close all
T1 = mri_resliced.transform;
T2 = mri_realigned.transform;

sourcemodel1 = ft_transform_geometry((T2/T1), sourcemodel);
headmodel1 = ft_transform_geometry((T2/T1), headmodel);

figure; hold on;
ft_plot_vol(headmodel1, 'facecolor', 'none'); alpha 0.5;
ft_plot_mesh(sourcemodel1, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
hold on;
ft_plot_headshape(headshape);

%%
sourcemodel = sourcemodel1;
headmodel = headmodel1;

%%
figure; hold on;
ft_plot_vol(headmodel, 'facecolor', 'none'); alpha 0.5;
ft_plot_mesh(sourcemodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
hold on;
ft_plot_headshape(headshape);

%%
% T1 = mri_resliced.transform;
% T2 = mri_realigned.transform;
% 
% % sourcemodel1 = ft_transform_geometry((T2/T1), sourcemodel);
% % headmodel1 = ft_transform_geometry((T2/T1), headmodel);
% 
% 
% pialdir = fullfile(indir,subj,'brainstorm_db/anat',subj,'tess_cortex_pial_low.mat');
% sourcemodel_bs = ft_read_headshape(pialdir);
% surface_sourcemodel = ft_convert_units(sourcemodel_bs, 'mm');
% 
% surface_sourcemodel = ft_transform_geometry((T2), surface_sourcemodel);

%%
% sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
% sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
% T   = tr_neuromag/tr_neuromagorg;
% sourcemodelT = ft_transform_geometry(T, surface_sourcemodel);
% InitTransf

% tr_neuromag = mri_realigned.transform;
% tr_ctf = mri_realigned_ctf.transform;
% tr_neuromagorg = mri_realigned.transformorig;


% figure; hold on
% ft_plot_vol(headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(surface_sourcemodel, 'edgecolor', 'k'); camlight
% view([90,0]);

% figure; hold on
% ft_plot_vol(headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(sourcemodel, 'edgecolor', 'k'); camlight
% view([90,0]);

%%

% ft_plot_ortho(mri_realigned.anatomy, 'transform', mri_realigned.transform, 'style', 'intersect');
% ft_plot_mesh(sourcemodelT, 'edgecolor', 'k'); camlight
% hold on
% ft_plot_vol(individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');

%%
% sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
% sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_ctf, surface_sourcemodel);
% sourcemodelT = ft_transform_geometry(tr_ctf, surface_sourcemodel);


%%
% transform = inv(ChannelMat.TransfMeg{strcmp('neuromag_head=>scs',ChannelMat.TransfMegLabels)});
%
% transform = [
%     0 -1 0 0;
%     1 0 0 0;
%     0 0 1 0;
%     0 0 0 1
%     ];
% 
% sourcemodelT = ft_transform_geometry(transform, sourcemodelT);

%%
% figure; hold on
% ft_plot_vol(headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% % ft_plot_mesh(sourcemodel_update, 'edgecolor', 'k'); camlight
% ft_plot_mesh(sourcemodel, 'edgecolor', 'k'); camlight
% view([90,0])
% % view([-90,0])
% % view([0,0])
% 
% %%
% % datain = t_data.all;
% datain = t_data.pst;


%%
% cfg = [];
% cfg.grad                = datain.grad;              % sensor positions
% % cfg.channel             = 'meggrad';                  % the used channels
% cfg.senstype            = 'meg';
% cfg.grid.pos            = sourcemodel.pos;           % source points
% cfg.grid.inside         = 1:size(sourcemodel.pos,1); % all source points are inside of the brain
% cfg.headmodel           = headmodel;              % volume conduction model
% 
% leadfield_mne = ft_prepare_leadfield(cfg,datain);
% 
% %%
% cfg                     = [];
% cfg.method              = 'mne';
% cfg.channel             = 'meggrad';
% cfg.senstype            = 'meg';
% cfg.grid                = leadfield_mne;
% cfg.headmodel           = headmodel;
% cfg.mne.prewhiten       = 'yes';
% cfg.mne.lambda          = 3;
% cfg.mne.scalesourcecov  = 'yes';
% source  = ft_sourceanalysis(cfg, datain);
% 
% %%
% % m = source.avg.pow;
% % [~,~,stats] = anova1(m, [],'off');
% 
% %%
% % data_cmb = ft_combineplanar([],datain); %Combine gradiometers
% % 
% % cfg = [];
% % cfg.layout = 'neuromag306cmb.lay';
% % figure; ft_multiplotER(cfg, data_cmb);
% 
% %%
% % peaksel = 3;
% % datain = t_data.pst;
% %%
% stat = stats.means;
% stat = (stat - min(stat(:))) ./ (max(stat(:)) - min(stat(:))); %
% %
% [psor,lsor] = findpeaks(stat,datain.time,'SortStr','descend');
% figure,plot(datain.time,stat),hold on
% text(lsor(1:peaksel),psor(1:peaksel),num2str((1:peaksel)'));
% grid
% box off
% xlabel('Time (sec)'); ylabel('Stats (normal)')
% set(gcf, 'Position',  [500, 500, 1200, 500]);
% 
% [~,lsor1] = findpeaks(stat,1:length(datain.time),'SortStr','descend');
% 
% savefig = fullfile(savepath,['sourceanova_',subj]);
% hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

%%
% source.tri = sourcemodelT.tri; %ft_sourceplot need the triangulation information, add to source reconstruction.
% 
% views =[180,0;0,0;90,0;180,90];
% for peaknum=1:peaksel
% %     figure,
% %     for i=1:peaksel
% %         subplot(2,2,i)
%         cfg = [];
%         cfg.method          = 'surface';
%         cfg.funparameter    = 'pow';
% %         cfg.funcolormap     = 'hot';
%         cfg.latency         = lsor(peaknum);     % The time-point to plot
%         cfg.colorbar        = 'no';
%         cfg.funcolormap        = brewermap(256, '*RdYlBu');
%         % cfg.avgovertime     = 'yes';
%         cfg.title           = [num2str(peaknum),': ',num2str(lsor(peaknum)),' Sec'];
%         ft_sourceplot(cfg, source);
% %         view([-100,20])
% %     savefig = fullfile(savepath,['MNE_peak',num2str(peaknum)]);
% %     hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
% end
% 
% %%
% data_cmb = ft_combineplanar([],t_data.pst); %Combine gradiometers
% 
% cfg = [];
% cfg.layout = 'neuromag306cmb.lay';
% figure; ft_multiplotER(cfg, data_cmb);
% 
% 
% %%
% source.tri = sourcemodelT.tri; %ft_sourceplot need the triangulation information, add to source reconstruction.
% 
% cfg = [];
% cfg.method          = 'surface';
% cfg.funparameter    = 'pow';
% cfg.funcolormap     = 'jet';
% cfg.latency         = [0.7 1];     % The time-point to plot
% cfg.colorbar        = 'no';
% cfg.avgovertime     = 'yes';
% ft_sourceplot(cfg, source);
% view([-100,20]);