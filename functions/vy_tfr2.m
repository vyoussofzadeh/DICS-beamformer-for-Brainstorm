function tfr = vy_tfr2(cfg_main, data)

data1 = data.pst;
% do tfr-decomposition
cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
% cfg.taper      = 'dpss';
cfg.foi        = .1:40;
cfg.keeptrials = 'yes';
cfg.t_ftimwin  = 3./cfg.foi;
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.tapsmofrq  = 0.8 *cfg.foi;
cfg.toi        = data1.time(1):0.05:data1.time(end);
tfr1        = ft_freqanalysis(cfg, data1);
if isempty(cfg_main.savefile) == 0
    save(cfg_main.savefile, 'tfr', '-v7.3');
end

data2 = data.bsl;
% do tfr-decomposition
cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
% cfg.taper      = 'dpss';
cfg.foi        = .1:40;
cfg.keeptrials = 'yes';
cfg.t_ftimwin  = 3./cfg.foi;
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.tapsmofrq  = 0.8 *cfg.foi;
cfg.toi        = data1.time(1):0.05:data1.time(end);
tfr2        = ft_freqanalysis(cfg, data1);
if isempty(cfg_main.savefile) == 0
    save(cfg_main.savefile, 'tfr', '-v7.3');
end

%%
tfr3 = tfr2;
tfr3.powspctrm = (tfr1.powspctrm / tfr2.powspctrm);
% ./(tfr1.powspctrm + tfr2.powspctrm);

%% Plotting TFR
cfg = [];
% cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'absolute';
% % cfg.showlabels   = 'yes';
% % cfg.funcolormap = flipud(brewermap(64,'RdBu'));
% cfg.zlim         = [-3e-27 3e-27];
cfg.layout       = cfg_main.lay;
figure
ft_multiplotTFR(cfg, tfr3);
% colorbar jet
% colormap(flipud(brewermap(64,'RdBu')));
if isempty(cfg_main.savefile) == 0
    hcp_write_figure([cfg_main.savefile,'.png'], gcf, 'resolution', 300);
end
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
