%% Volumetric-based analysis
mridir = fullfile(indir,subj,'brainstorm_db/anat');
d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
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
pialdir = fullfile(indir,subj,'brainstorm_db/anat',subj,'tess_cortex_pial_low.mat');
sourcemodel = ft_read_headshape(pialdir);
surface_sourcemodel = ft_convert_units(sourcemodel, 'mm');

%%
% sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
% sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
% T   = tr_neuromag/tr_neuromagorg;
% sourcemodelT = ft_transform_geometry(T, surface_sourcemodel);
% InitTransf

% tr_neuromag = mri_realigned.transform;
tr_ctf = mri_realigned_ctf.transform;
tr_neuromagorg = mri_realigned.transformorig;

% 
% figure; hold on
% ft_plot_vol(individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(sourcemodelT, 'edgecolor', 'k'); camlight
% view([90,0]);

%%

% ft_plot_ortho(mri_realigned.anatomy, 'transform', mri_realigned.transform, 'style', 'intersect');
% ft_plot_mesh(sourcemodelT, 'edgecolor', 'k'); camlight
% hold on
% ft_plot_vol(individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');

%%
% sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_ctf, surface_sourcemodel);
% sourcemodelT = ft_transform_geometry(tr_ctf, surface_sourcemodel);


%%
% transform = inv(ChannelMat.TransfMeg{strcmp('neuromag_head=>scs',ChannelMat.TransfMegLabels)});
%
transform = [
    0 -1 0 0;
    1 0 0 0;
    0 0 1 0;
    0 0 0 1
    ];

sourcemodelT = ft_transform_geometry(transform, sourcemodelT);

%%
figure; hold on
ft_plot_vol(outanat.individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(sourcemodel_update, 'edgecolor', 'k'); camlight
ft_plot_mesh(surface_sourcemodel, 'edgecolor', 'k'); camlight
view([90,0])
% view([-90,0])
% view([0,0])

%%
% datain = t_data.all;
datain = t_data.pst;


%%
cfg = [];
cfg.grad                = datain.grad;              % sensor positions
% cfg.channel             = 'meggrad';                  % the used channels
cfg.senstype            = 'meg';
cfg.grid.pos            = sourcemodelT.pos;           % source points
cfg.grid.inside         = 1:size(sourcemodelT.pos,1); % all source points are inside of the brain
cfg.headmodel           = individual_headmodel;              % volume conduction model

leadfield_mne = ft_prepare_leadfield(cfg,datain);

%%
cfg                     = [];
cfg.method              = 'mne';
cfg.channel             = 'meggrad';
cfg.senstype            = 'meg';
cfg.grid                = leadfield_mne;
cfg.headmodel           = individual_headmodel;
cfg.mne.prewhiten       = 'yes';
cfg.mne.lambda          = 3;
cfg.mne.scalesourcecov  = 'yes';
source  = ft_sourceanalysis(cfg, datain);

%%
m = source.avg.pow;
[~,~,stats] = anova1(m, [],'off');

%%
% data_cmb = ft_combineplanar([],datain); %Combine gradiometers
% 
% cfg = [];
% cfg.layout = 'neuromag306cmb.lay';
% figure; ft_multiplotER(cfg, data_cmb);

%%
peaksel = 3;
% datain = t_data.pst;
%%
stat = stats.means;
stat = (stat - min(stat(:))) ./ (max(stat(:)) - min(stat(:))); %
%
[psor,lsor] = findpeaks(stat,datain.time,'SortStr','descend');
figure,plot(datain.time,stat),hold on
text(lsor(1:peaksel),psor(1:peaksel),num2str((1:peaksel)'));
grid
box off
xlabel('Time (sec)'); ylabel('Stats (normal)')
set(gcf, 'Position',  [500, 500, 1200, 500]);

[~,lsor1] = findpeaks(stat,1:length(datain.time),'SortStr','descend');

savefig = fullfile(savepath,['sourceanova_',subj]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

%%
source.tri = sourcemodelT.tri; %ft_sourceplot need the triangulation information, add to source reconstruction.

views =[180,0;0,0;90,0;180,90];
for peaknum=1:peaksel
%     figure,
%     for i=1:peaksel
%         subplot(2,2,i)
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
%         cfg.funcolormap     = 'hot';
        cfg.latency         = lsor(peaknum);     % The time-point to plot
        cfg.colorbar        = 'no';
        cfg.funcolormap        = brewermap(256, '*RdYlBu');
        % cfg.avgovertime     = 'yes';
        cfg.title           = [num2str(peaknum),': ',num2str(lsor(peaknum)),' Sec'];
        ft_sourceplot(cfg, source);
%         view([-100,20])
%     savefig = fullfile(savepath,['MNE_peak',num2str(peaknum)]);
%     hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
end

%%
data_cmb = ft_combineplanar([],t_data.pst); %Combine gradiometers

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
figure; ft_multiplotER(cfg, data_cmb);


%%
source.tri = sourcemodelT.tri; %ft_sourceplot need the triangulation information, add to source reconstruction.

cfg = [];
cfg.method          = 'surface';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';
cfg.latency         = [0.7 1];     % The time-point to plot
cfg.colorbar        = 'no';
cfg.avgovertime     = 'yes';
ft_sourceplot(cfg, source);
view([-100,20]);