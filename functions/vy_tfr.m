function tfr = vy_tfr(cfg_main, data)

% do tfr-decomposition
cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
% cfg.taper      = 'dpss';
cfg.foi        = 1:3:40;
cfg.keeptrials = 'yes';
cfg.t_ftimwin  = 3./cfg.foi;
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.tapsmofrq  = 0.8 *cfg.foi;
cfg.toi        = -0.5:0.05:2;
tfr        = ft_freqanalysis(cfg, data);
set(gcf,'name',cfg_main.subj,'numbertitle','off')
if isempty(cfg_main.savefile) == 0
    save(cfg_main.savefile, 'tfr', '-v7.3');
end

%% Plotting TFR
% cfg = [];
% cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'absolute';
% % cfg.showlabels   = 'yes';
% % cfg.funcolormap = flipud(brewermap(64,'RdBu'));
% cfg.zlim         = [-3e-27 3e-27];
% cfg.layout       = cfg_main.lay;
% figure
% ft_multiplotTFR(cfg, tfr);
% set(gcf,'name',cfg_main.subj,'numbertitle','off')
% % colorbar jet
% % colormap(flipud(brewermap(64,'RdBu')));
% if isempty(cfg_main.savefile) == 0
%     hcp_write_figure([cfg_main.savefile,'.png'], gcf, 'resolution', 300);
% end

%%
% figure
% cfg                = [];
% cfg.baseline       = [-.5 -0.1]; % a 3 s baseline around the event as it has no clear start or end.
% cfg.baselinetype   = 'absolute';
% cfg.zlim         = [-3e-27 3e-27];
% cfg.layout       = lay;
% ft_singleplotTFR(cfg,tfr);
% % colormap(flipud(brewermap(64,'RdBu')));
% 
% if isempty(savepath) == 0
%     hcp_write_figure([savepath,'_2.png'], gcf, 'resolution', 300);
% end
