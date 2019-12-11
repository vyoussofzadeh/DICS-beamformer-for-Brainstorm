

mtag = 'dics'; outd.vol = fullfile(outd.sub,'source_surface',mtag);
cfg = [];
cfg.grid = leadfield_mne;
cfg.headmodel = individual_headmodel;
cfg.sens = sens;
cfg.outputdir = outd.vol;
% cfg.template_grid = template_grid;
% cfg.template_mri = template_mri;
source_diff_dics = vy_source_dics_surface(cfg,ep_data);

%%
source_diff_dics.tri = sourcemodelT.tri; %ft_sourceplot need the triangulation information, add to source reconstruction.


% cfg = [];
% cfg.method          = 'ortho';
% cfg.funparameter    = 'pow';
% cfg.funcolormap     = 'jet';
% cfg.latency         = 0;     % The time-point to plot
% cfg.colorbar        = 'no';
% % cfg.avgovertime     = 'yes';
% ft_sourceplot(cfg, source);
% view([-100,20]);

%%
% source_diff_dics.time = 1;
% % 
% % source_diff_dics.avg.pow = source_diff_dics.pow;
% cfg = [];
% cfg.funparameter = 'pow';
% cfg.method = 'ortho';  % plot slices
% ft_sourceplot(cfg, source);


%%
% views =[180,0;0,0;90,0;180,90];
% 
% bnd = [];
% bnd.pnt = sourcemodelT.pos;
% bnd.tri = sourcemodelT.tri;
% bnd.funcolormap =  brewermap(256, '*RdYlBu');
% m1 = source_diff_dics.pow;
% m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); %
% figure,
% ft_plot_mesh(bnd, 'vertexcolor', m1, 'maskstyle', 'opacity');
% colormap(brewermap(256, '*RdYlBu'));
% view(views(i,:))


%%
% views =[-90,0; 90,0; 180, 0; 0,90];
% bnd = [];
% bnd.pnt = sourcemodelT.pos;
% bnd.tri = sourcemodelT.tri;
% m1 = source_diff_dics.pow;
% % m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); %
% 
% figure,
% for i=1:4
%     subplot(2,2,i)
%     ft_plot_mesh(bnd, 'vertexcolor', m1, 'maskstyle', 'opacity');
%     view(views(i,:))
% end
% colorbar
% colormap(brewermap(256, '*RdYlBu'));
% mtit([num2str(peaknum),': ',num2str(lsor(peaknum)),' Sec'],'fontsize',14,'color',[0 0 0],'xoff',0,'yoff',0);
% savefig = fullfile(savepath,['MNE_peak',num2str(peaknum)]);
% hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

%%

source_diff_dics.time = 1;
cfg = [];
cfg.method          = 'surface';
% cfg.method        = 'slice';
cfg.funparameter    = 'pow';
% cfg.maskparameter  = cfg.funparameter;
m1 = source_diff_dics.pow;
m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); 
source_diff_dics.pow = m1;
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
ft_sourceplot(cfg, source_diff_dics);
colorbar
view([-100,0])
light ('Position',[-70 20 50])
material dull

%%
% cfg = [];
% cfg.projectmom = 'yes';
% sdFC  = ft_sourcedescriptives(cfg,source_diff_dics);
% 
% cfg = [];
% cfg.funparameter = 'pow';
% ft_sourceplot(cfg,sdFC);

