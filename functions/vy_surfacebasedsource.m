function  [source, individual_headmodel, individual_grid, sourcemodel] = vy_surfacebasedsource(cfg_main)

if exist(fullfile(cfg_main.outputmridir,['surfanat_',cfg_main.subj, '_', cfg_main.task, '.mat']), 'file') == 2
    load(fullfile(cfg_main.outputmridir,['surfanat_',cfg_main.subj, '_', cfg_main.task, '.mat']));
else
    
    % BS2FT
    pialdir = fullfile(cfg_main.datadir,cfg_main.subj,'brainstorm_db/anat',cfg_main.subj,'tess_cortex_pial_low.mat');
    sourcemodel = ft_read_headshape(pialdir);
    
    %     figure
    %     ft_plot_mesh(sourcemodel, 'vertexcolor', sourcemodel.curv)
    
    [individual_headmodel, individual_grid] = vy_bs2ft_headmodel(cfg_main.data, sourcemodel);
    
    figure;
    ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    hold on;
    ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
    view ([-10 40 10])
    
    % hsfile = datafile; % headshape
    % headshape = ft_read_headshape(hsfile);
    % headshape = ft_convert_units(headshape,'mm');
    % ft_plot_headshape(headshape);
    save(fullfile(cfg_main.outputmridir,['surfanat_',cfg_main.subj, '_', cfg_main.task, '.mat']), 'individual_grid','individual_headmodel','sourcemodel');
end

%%
% figure; ft_plot_mesh(individual_grid, 'edgecolor', 'k'); camlight
% 
% figure; hold on
% ft_plot_vol(individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(individual_grid, 'edgecolor', 'k'); camlight


%%
savepath = fullfile(cfg_main.outputdir);
if exist(savepath, 'file') == 0, mkdir(savepath), end

%% MNE
cfg               = [];
cfg.method        = 'mne';
cfg.grid          = individual_grid;
cfg.headmodel     = individual_headmodel;
cfg.mne.lambda    = 3;
cfg.mne.scalesourcecov = 'yes';
% cfg.mne.normalize = 'yes';
% cfg.mne.prewhiten = 'yes';
% cfg.mne.keepfilter= 'yes';
% cfg.mne.projectnoise = 'yes';
% cfg.projectmom = 'yes';
source            = ft_sourceanalysis(cfg,cfg_main.data);
m = source.avg.pow;
[~,~,stats] = anova1(m, [],'off');

% figure, plot(cfg_main.data.time,stats.means), xlabel('time'); ylabel('Stats')

%%
% remove center of head bias by using the neural activity index (NAI)
% save into single structure to save memory
source = ft_struct2single(source);

%%
% cfg = [];
% cfg.projectmom = 'yes';
% msource  = ft_sourcedescriptives(cfg,source);

%%
% cfg = [];
% cfg.funparameter = 'pow';
% ft_sourcemovie(cfg,source);

%%
% cfg=[];
% cfg.method='eloreta';
% cfg.headmodel=individual_headmodel;
% cfg.grid=individual_grid;
% cfg.projectnoise='yes';
% cfg.keepcsd='yes';
% cfg.eloreta.projectnoise='yes';
% cfg.eloreta.keepcsd='yes';
% cfg.eloreta.keepmom='yes';
% cfg.eloreta.lambda=1e-5;
% cfg.lambda=cfg.eloreta.lambda;
% source=ft_sourceanalysis(cfg,cfg_main.data);

%%
peaksel = cfg_main.peaksel;

%%
stat = stats.means;
stat = (stat - min(stat(:))) ./ (max(stat(:)) - min(stat(:))); %
%
[psor,lsor] = findpeaks(stat,cfg_main.data.time,'SortStr','descend');
figure,plot(cfg_main.data.time,stat),hold on
text(lsor(1:peaksel),psor(1:peaksel),num2str((1:peaksel)'));
grid
box off
xlabel('Time (sec)'); ylabel('Stats (normal)')
set(gcf, 'Position',  [500, 500, 1200, 500]);

[~,lsor1] = findpeaks(stat,1:length(cfg_main.data.time),'SortStr','descend');

savefig = fullfile(savepath,['sourceanova_',cfg_main.subj]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

%%
views =[180,0;0,0;90,0;180,90];
for peaknum=1:peaksel
    figure,
    for i=1:peaksel
        subplot(2,2,i)
        bnd.pnt = individual_grid.pos;
        bnd.tri = individual_grid.tri;
        bnd.funcolormap =  brewermap(256, '*RdYlBu');
        m1 = m(:,lsor1(peaknum));
        m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); %
        ft_plot_mesh(bnd, 'vertexcolor', m1, 'maskstyle', 'opacity');
        view(views(i,:))
    end
    % colorbar
    colormap(brewermap(256, '*RdYlBu'));
    mtit([num2str(peaknum),': ',num2str(lsor(peaknum)),' Sec'],'fontsize',14,'color',[0 0 0],'xoff',0,'yoff',0);
    savefig = fullfile(savepath,['MNE_peak',num2str(peaknum)]);
    hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
end


%%
% figure,
% for i=1:3
%     subplot(2,2,i)
%     % m1 = mean(m,2); % plotting the result at 400 ms
%     bnd.pnt = individual_grid.pos;
%     bnd.tri = individual_grid.tri;
%     % bnd.funcolorlim   = [0.1 1];
%     % bnd.opacitylim    = [0.1 1];
%     bnd.funcolormap =  brewermap(256, '*RdYlBu');
%     m1 = source.avg.pow(:,lsor(i)); % plotting the result at 400 ms
%     m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); %
%     ft_plot_mesh(bnd, 'vertexcolor', m1, 'maskstyle', 'opacity');
%         colorbar
%         colormap(brewermap(256, '*RdYlBu')); title(['peak: ', num2str(i)])
% end

%%
% [mx, idx] = max(stats.means);
%
% m1 = source.avg.pow(:,idx); % plotting the result at 400 ms
% % m1 = source.avg.pow(:,651); % plotting the result at 400 ms
% % m1 = source.avg.pow(:,110); % plotting the result at 400 ms
% % m1 = source.avg.pow(:,910); % plotting the result at 400 ms
% % m1 = source.avg.pow(:,545); % plotting the result at 400 ms
%
% m1 = (m1 - min(m1(:))) ./ (max(m1(:)) - min(m1(:))); %
%
% figure
% % m1 = mean(m,2); % plotting the result at 400 ms
% bnd.pnt = individual_grid.pos;
% bnd.tri = individual_grid.tri;
% % bnd.funcolorlim   = [0.1 1];
% % bnd.opacitylim    = [0.1 1];
% bnd.funcolormap =  brewermap(256, '*RdYlBu');
% ft_plot_mesh(bnd, 'vertexcolor', m1, 'maskstyle', 'opacity');
% colorbar
% colormap(brewermap(256, '*RdYlBu'));


%% BF
% method = 'lcmv';
% cfg                  = [];
% cfg.method           = 'lcmv';
% cfg.grid             = individual_grid; % leadfield, which has the grid information
% cfg.headmodel        = individual_headmodel; % volume conduction model (headmodel)
% cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
% cfg.lcmv.lambda      = '10%';
% cfg.lcmv.keepfilter  = 'yes';
% sourceAll = ft_sourceanalysis(cfg, t_data.app);
% cfg.grid.filter = sourceAll.avg.filter;
% s_data.bsl = ft_sourceanalysis(cfg, t_data.bsl);
% s_data.pst = ft_sourceanalysis(cfg, t_data.pst);
%
% s_data2.bsl  = ft_sourcedescriptives([], s_data.bsl); % to get the neural-activity-index
% s_data2.pst  = ft_sourcedescriptives([], s_data.pst); % to get the neural-activity-index
%
% %
% cfg = [];
% cfg.parameter = 'pow';
% %     cfg.operation = '(x1-x2)/(x1+x2)';
% cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
% source_ratio = ft_math(cfg,s_data.pst,s_data.bsl);
% %     source_ratio.pow(source_ratio.pow>0)=0;
% %     source_ratio.pow = abs(source_ratio.pow);
%
% %
% figure
% m = source_ratio.pow;
% bnd.pnt = individual_grid.pos;
% bnd.tri = individual_grid.tri;
% ft_plot_mesh(bnd, 'vertexcolor', m);
% colorbar


