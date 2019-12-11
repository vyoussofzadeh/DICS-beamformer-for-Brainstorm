% d = rdir(fullfile(bsdatadir,subj,['/@*',tag,'*'],'/channel_vectorview306*.mat'));
% if isempty(d)
%     d = rdir(fullfile(bsdatadir,subj,['/*',tag,'*'],'/channel_vectorview306*.mat'));
%     if exist('tag1','var') && isempty(d)
%         d = rdir(fullfile(bsdatadir,subj,['/*',tag1,'*'],'/channel_vectorview306*.mat'));
%     end
% end
% 
% %%
% clear datafilebf
% for i=1:length(d)
%     [pathstr, name] = fileparts(d(i).name);
%     %     databf{i} = pathstr;
%     datafilebf{i} = d(i).name;
% end
% disp(datafilebf')
% channel = load(datafilebf{1});
% disp('=============')
% disp(datafilebf{1})
% disp('was loaded')

%%
bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
bsanatdir = fullfile(indir,subj,'brainstorm_db/anat');

%%
d = rdir(fullfile(bsdatadir,subj,['/Run*',tag,'*'],'headmodel_surf_os_meg.mat'));

clear HeadModelFile
for i=1:length(d) 
    [pathstr, name] = fileparts(d(i).name); HeadModelFile{i} = d(i).name;
end
disp(HeadModelFile')

% fp = '/MEG_data/epilepsy/hulan_eleanor/brainstorm_db/data/hulan_eleanor/DFN_hulan_eleanor_IC_data';
% HeadModelFile = fullfile(fp,'/headmodel_surf_os_meg.mat');

SurfaceFile = fullfile(bsanatdir,subj,'tess_cortex_pial_low.mat');

% SurfaceFile = fullfile(indir,subj,'brainstorm_db/anat',subj,'tess_cortex_pial_low.mat');
sourcemodel = ft_read_headshape(SurfaceFile);
headshape = ft_read_headshape(datafile);headshape = ft_convert_units(headshape,'m');

% close all
% set(0,'DefaultFigureWindowStyle','normal')
% addpath('/usr/local/MATLAB_Tools/brainstorm3')
% brainstorm
% HeadModelMat = in_bst_headmodel(HeadModelFile{1});

d = rdir(fullfile(bsdatadir,subj,['/Run*',tag,'*'],'channel_vectorview306_acc1.mat'));

clear HeadModelFile AllChannelFiles
for i=1:length(d) 
    [pathstr, name] = fileparts(d(i).name); 
    AllChannelFiles{i} = d(i).name;
end
disp(AllChannelFiles')
% 
% AllChannelFiles = fullfile(fp,'channel_vectorview306_acc1.mat');
ChannelMat = load(AllChannelFiles{1});
% ChannelMat = in_bst_channel(AllChannelFiles{1});
% % 
% iChannelsData = 1:length(ChannelMat.Channel);
% % 
% [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);

%%

% seg = load(fullfile(bsanatdir,subj,'tess_aseg.mat'));
subjectimage_T1 = load(fullfile(bsanatdir,subj,'subjectimage_T1.mat'));

nii_mri_org = 'MRI_org.nii';
cfg                 = [];
cfg.filename        = nii_mri_org;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, subjectimage_T1);


cfg = [];
cfg.output = 'brain';
seg = ft_volumesegment(cfg, subjectimage_T1);

cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 10000;
bnd = ft_prepare_mesh(cfg, individual_seg);

cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, bnd);

%% BS
% brainstorm
% HeadModelMat = in_bst_headmodel(HeadModelFile{1});
% 
% ChannelMat = in_bst_channel(AllChannelFiles{1});
% % 
% iChannelsData = 1:length(ChannelMat.Channel);
% % 
[ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
% %
% brainstorm stop
% 
% figure; hold on;
% ft_plot_vol(ftHeadmodel); alpha 0.5;
% ft_plot_mesh(ftLeadfield, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
% % ft_plot_mesh(leadfield, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
% hold on;
% ft_plot_headshape(sourcemodel);
% hold on;
% ft_plot_headshape(headshape1);
%

%%
% figure; hold on;
% % ft_plot_vol(ftHeadmodel); alpha 0.5;
% ft_plot_mesh(ftLeadfield, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
% hold on;
% ft_plot_headshape(sourcemodel);
% a = ChannelMat.HeadPoints;
% a = [];
% a.pos = ChannelMat.HeadPoints.Loc';
% a.label = ChannelMat.HeadPoints.Label;
% ft_plot_headshape(a);

headshape1 = [];
headshape1 = ft_transform_geometry(ChannelMat.TransfMeg{2},headshape);
% 
% figure; hold on;
% % ft_plot_vol(ftHeadmodel); alpha 0.5;
% ft_plot_mesh(ftLeadfield, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
% hold on;
% ft_plot_headshape(sourcemodel);
% ft_plot_headshape(headshape1);

%%
cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, sourcemodel);

%%
% figure; hold on;
% ft_plot_vol(headmodel); alpha 0.5;
% % ft_plot_mesh(ftLeadfield, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
% hold on;
% ft_plot_headshape(sourcemodel);

%%
ep_data1 = ep_data;
ep_data1.all.grad = ft_transform_geometry(ChannelMat.TransfMeg{2},ep_data1.all.grad);
ep_data1.app.grad = ft_transform_geometry(ChannelMat.TransfMeg{2},ep_data1.app.grad);
ep_data1.bsl.grad = ft_transform_geometry(ChannelMat.TransfMeg{2},ep_data1.bsl.grad);
ep_data1.pst.grad = ft_transform_geometry(ChannelMat.TransfMeg{2},ep_data1.pst.grad);

%%
datain = t_data.all;
cfg = [];
cfg.grad                = ep_data1.all.grad;              % sensor positions
% cfg.channel             = 'meggrad';                  % the used channels
cfg.senstype            = 'meg';
% cfg.sourcemodel = sourcemodel;
cfg.grid.pos            = sourcemodel.pos;           % source points
cfg.grid.inside         = 1:size(sourcemodel.pos,1); % all source points are inside of the brain
cfg.headmodel           = headmodel;              % volume conduction model
leadfield = ft_prepare_leadfield(cfg,datain);


% create leadfield
% cfg                  = [];
% cfg.grad             = ep_data1.all.grad;  % gradiometer distances
% cfg.headmodel        = headmodel;   % volume conduction headmodel
% cfg.sourcemodel      = sourcemodel;
% cfg.channel          = {'meg'};
% cfg.singleshell.batchsize = 2000;
% leadfield            = ft_prepare_leadfield(cfg);

%%
figure; hold on;
ft_plot_vol(headmodel); 
alpha 0.5;
ft_plot_mesh(leadfield, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
hold on;
ft_plot_headshape(sourcemodel);
hold on;
ft_plot_headshape(headshape1);

%%
figure; hold on;
ft_plot_vol(headmodel); alpha 0.5;
ft_plot_mesh(leadfield, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
hold on;
ft_plot_headshape(sourcemodel);
hold on;
ft_plot_headshape(headshape1);
% ft_plot_sens(ep_data1.pst.grad, 'unit', 'mm', 'coilsize', 10, 'chantype', 'meggrad');
grad = ft_convert_units(ep_data1.pst.grad,'m');ft_plot_sens(grad)

%%
% figure,
% ft_plot_mesh(headmodel.bnd);

%%
sens = ft_read_sens(datafile,'senstype','meg');
sens = ft_convert_units(sens,'m');

sens1 = ft_transform_geometry(ChannelMat.TransfMeg{2}, sens);

figure,
ft_plot_sens(sens1)
hold on
ft_plot_headshape(sourcemodel);

%%
cfg = [];
cfg.headmodel = headmodel;
cfg.sourcemodel = sourcemodel;
cfg.leadfield   = leadfield;
cfg.mtag = 'dics_fs';
cfg.sens = sens;
cfg.subj = subj;
cfg.outputdir = outd.vol;
vy_source_dics_fs(cfg, ep_data);

%%
cfg = [];
cfg.grid = leadfield;
cfg.headmodel = headmodel;
cfg.sens = sens1;
cfg.outputdir = [];
cfg.sourcmodel = sourcemodel;
% cfg.template_grid = template_grid;
% cfg.template_mri = template_mri;
source_diff_dics = vy_source_dics_surface(cfg,ep_data1);

%%
figure;
bnd.pnt = sourcemodel.pos;
bnd.tri = sourcemodel.tri;
bnd.funcolormap =  brewermap(256, '*RdYlBu');
cfg.opacitymap    = 'rampdown';
ft_plot_mesh(bnd, 'vertexcolor', source_diff_dics.pow);
colorbar
light ('Position',[-70 20 50])
material dull


%%
figure; hold on
ft_plot_vol(headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(sourcemodel, 'edgecolor', 'k'); camlight

figure; ft_plot_mesh(sourcemodel, 'edgecolor', 'k'); camlight 


%%
% source_diff_dics.tri = sourcemodel.tri;
% source_diff_dics.tri = headmodel.bnd.tri;

source_diff_dics1 = source_diff_dics;
source_diff_dics1.time = 1;

cfg = [];
cfg.method          = 'surface';
% cfg.method        = 'slice';
cfg.funparameter    = 'pow';
% cfg.maskparameter  = cfg.funparameter;
m1 = source_diff_dics1.pow;
m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); 
source_diff_dics1.pow = m1;
%         cfg.funcolormap     = 'hot';
cfg.latency         = 1;     % The time-point to plot
cfg.colorbar        = 'no';
cfg.funcolormap     = brewermap(256, '*RdYlBu');
% cfg.funcolorlim  = [-1, 0];
% cfg.funcolorlim   = [0 0.5];
% cfg.opacitylim    = [0 0.5];

% cfg.funcolorlim    = [0.0 1.2];
% cfg.funcolormap    = 'jet';
% cfg.opacitylim     = [0.0 1.2];
cfg.opacitymap    = 'rampdown';
% cfg.surffile       = 'surface_white_both.mat';
% cfg.projmethod     = 'nearest'; 
% cfg.surffile       = 'surface_white_both.mat';
% cfg.surfdownsample = 10;
% cfg.surffile = 'surface_white_both.mat';
% cfg.surfinflated = 'surface_inflated_both_caret.mat';
% cfg.surfdownsample = 10;
% cfg.avgovertime     = 'yes';
% cfg.title           = [num2str(peaknum),': ',num2str(lsor(peaknum)),' Sec'];
% cfg.projmethod     = 'nearest';
% cfg.surffile   = 'surface_inflated_both_caret.mat';
% cfg.projthresh     = 0.8;
ft_sourceplot(cfg, source_diff_dics1);
colorbar
% view([-100,0])
light ('Position',[-70 20 50])
material dull

%%
views =[180,0;0,0;90,0;180,90];
for peaknum=1:1
    figure,
    for i=1:1
        %         subplot(2,2,i)
        bnd.pnt = sourcemodel.pos;
        bnd.tri = sourcemodel.tri;
        bnd.funcolormap =  brewermap(256, '*RdYlBu');
        %         m1 = m(:,lsor1(peaknum));
        %         m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); %
        %         m1 =
        ft_plot_mesh(bnd, 'vertexcolor', m1, 'maskstyle', 'opacity');
        %         view(views(i,:))
    end
    % colorbar
    colormap(brewermap(256, '*RdYlBu'));
    %     mtit([num2str(peaknum),': ',num2str(lsor(peaknum)),' Sec'],'fontsize',14,'color',[0 0 0],'xoff',0,'yoff',0);
    %     savefig = fullfile(savepath,['MNE_peak',num2str(peaknum)]);
    %     hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
end

%%
DataFile = fullfile(fp,'data_DFN_hulan_eleanor_IC_average_190806_1107.mat');
% DataFile = sInputs(iChanInputs(iInput)).FileName;
DataMat = in_bst_data(DataFile);
ftData = out_fieldtrip_data(DataMat, ChannelMat, iChannelsData, 1);

cfg           = [];
%             cfg.grid      = ftLeadfield;
cfg.sourcemodel      = ftLeadfield;
cfg.headmodel = ftHeadmodel;
cfg.method      = 'lcmv';
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.lambda      = '10%';
cfg.lcmv.keepfilter  = 'yes';
ftSource = ft_sourceanalysis(cfg, ftData);

% source = ft_struct2single(ftSource);

pow = sqrt(ftSource.avg.pow);

%%
views =[180,0;0,0;90,0;180,90];
for peaknum=1:1
    figure,
    for i=1:1
        %         subplot(2,2,i)
        bnd.pnt = sourcemodel.pos;
        bnd.tri = sourcemodel.tri;
        bnd.funcolormap =  brewermap(256, '*RdYlBu');
        %         m1 = m(:,lsor1(peaknum));
        %         m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); %
        %         m1 =
        ft_plot_mesh(bnd, 'vertexcolor', m1, 'maskstyle', 'opacity');
        %         view(views(i,:))
    end
    % colorbar
    colormap(brewermap(256, '*RdYlBu'));
    %     mtit([num2str(peaknum),': ',num2str(lsor(peaknum)),' Sec'],'fontsize',14,'color',[0 0 0],'xoff',0,'yoff',0);
    %     savefig = fullfile(savepath,['MNE_peak',num2str(peaknum)]);
    %     hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
end


%%
source_diff_dics1 = ftSource;
source_diff_dics1.tri  = sourcemodel.tri;
source_diff_dics1.time = 1;
cfg = [];
cfg.method          = 'surface';
% cfg.method        = 'slice';
cfg.funparameter    = 'pow';
% cfg.maskparameter  = cfg.funparameter;
m1 = source_diff_dics1.avg.pow;
m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); 
source_diff_dics1.avg.pow = m1;
%         cfg.funcolormap     = 'hot';
cfg.latency         = 1;     % The time-point to plot
cfg.colorbar        = 'no';
cfg.funcolormap     = brewermap(256, '*RdYlBu');
% cfg.funcolorlim  = [-1, 0];
% cfg.funcolorlim   = [0 0.5];
% cfg.opacitylim    = [0 0.5];

% cfg.funcolorlim    = [0.0 1.2];
% cfg.funcolormap    = 'jet';
% cfg.opacitylim     = [0.0 1.2];
cfg.opacitymap    = 'rampdown';
% cfg.surffile       = 'surface_white_both.mat';
% cfg.projmethod     = 'nearest'; 
% cfg.surffile       = 'surface_white_both.mat';
% cfg.surfdownsample = 10;
% cfg.surffile = 'surface_white_both.mat';
% cfg.surfinflated = 'surface_inflated_both_caret.mat';
% cfg.surfdownsample = 10;
% cfg.avgovertime     = 'yes';
% cfg.title           = [num2str(peaknum),': ',num2str(lsor(peaknum)),' Sec'];
% cfg.projmethod     = 'nearest';
% cfg.surffile   = 'surface_inflated_both_caret.mat';
% cfg.projthresh     = 0.8;
ft_sourceplot(cfg, source_diff_dics1);
colorbar
view([-100,0])
light ('Position',[-70 20 50])
material dull




% %% Volumetric-based analysis
% mridir = fullfile(indir,subj,'brainstorm_db/anat');
% d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
% clear fid
% if ~isempty(d)
%     sMRI1 = d.name;
%     load(sMRI1);
%     fid.SCS = SCS;
%     fid.NCS = NCS;
%     mripfile = fullfile(mridir,'T1.nii');
%     if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
%     
%     cfg = [];
%     cfg.megdata = t_data.app;
%     cfg.mripfile = mripfile;
%     cfg.hsfile = datafile; % headshape;
%     cfg.fid = fid;
%     cfg.outputmridir = outputmridir;
%     cfg.subj = subj;
%     cfg.plotflag = 2;
%     outanat = vy_mri_neuromag2(cfg);
%     
% end
% 
% %%
% pialdir = fullfile(indir,subj,'brainstorm_db/anat',subj,'tess_cortex_pial_low.mat');
% sourcemodel = ft_read_headshape(pialdir);
% surface_sourcemodel = ft_convert_units(sourcemodel, 'mm');
% 
% %%
% % sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
% % sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
% % T   = tr_neuromag/tr_neuromagorg;
% % sourcemodelT = ft_transform_geometry(T, surface_sourcemodel);
% % InitTransf
% 
% % tr_neuromag = mri_realigned.transform;
% tr_ctf = mri_realigned_ctf.transform;
% tr_neuromagorg = mri_realigned.transformorig;
% 
% % 
% % figure; hold on
% % ft_plot_vol(individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% % ft_plot_mesh(sourcemodelT, 'edgecolor', 'k'); camlight
% % view([90,0]);
% 
% %%
% 
% % ft_plot_ortho(mri_realigned.anatomy, 'transform', mri_realigned.transform, 'style', 'intersect');
% % ft_plot_mesh(sourcemodelT, 'edgecolor', 'k'); camlight
% % hold on
% % ft_plot_vol(individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% 
% %%
% % sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_neuromag, surface_sourcemodel);
% sourcemodelT = ft_transform_geometry(tr_neuromagorg/tr_ctf, surface_sourcemodel);
% % sourcemodelT = ft_transform_geometry(tr_ctf, surface_sourcemodel);
% 
% 
% %%
% % transform = inv(ChannelMat.TransfMeg{strcmp('neuromag_head=>scs',ChannelMat.TransfMegLabels)});
% %
% transform = [
%     0 -1 0 0;
%     1 0 0 0;
%     0 0 1 0;
%     0 0 0 1
%     ];
% 
% sourcemodelT = ft_transform_geometry(transform, sourcemodelT);
% 
% %%
% figure; hold on
% ft_plot_vol(outanat.individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% % ft_plot_mesh(sourcemodel_update, 'edgecolor', 'k'); camlight
% ft_plot_mesh(surface_sourcemodel, 'edgecolor', 'k'); camlight
% view([90,0])
% % view([-90,0])
% % view([0,0])
% 
% %%
% % datain = t_data.all;
% datain = t_data.pst;
% 
% 
% %%
% cfg = [];
% cfg.grad                = datain.grad;              % sensor positions
% % cfg.channel             = 'meggrad';                  % the used channels
% cfg.senstype            = 'meg';
% cfg.grid.pos            = sourcemodelT.pos;           % source points
% cfg.grid.inside         = 1:size(sourcemodelT.pos,1); % all source points are inside of the brain
% cfg.headmodel           = individual_headmodel;              % volume conduction model
% 
% leadfield_mne = ft_prepare_leadfield(cfg,datain);
% 
% %%
% cfg                     = [];
% cfg.method              = 'mne';
% cfg.channel             = 'meggrad';
% cfg.senstype            = 'meg';
% cfg.grid                = leadfield_mne;
% cfg.headmodel           = individual_headmodel;
% cfg.mne.prewhiten       = 'yes';
% cfg.mne.lambda          = 3;
% cfg.mne.scalesourcecov  = 'yes';
% source  = ft_sourceanalysis(cfg, datain);
% 
% %%
% m = source.avg.pow;
% [~,~,stats] = anova1(m, [],'off');
% 
% %%
% % data_cmb = ft_combineplanar([],datain); %Combine gradiometers
% % 
% % cfg = [];
% % cfg.layout = 'neuromag306cmb.lay';
% % figure; ft_multiplotER(cfg, data_cmb);
% 
% %%
% peaksel = 3;
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
% 
% %%
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