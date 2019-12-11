function comp = vy_pca(data,lay)

% cfg            = [];
% cfg.method     = 'runica';
% cfg.numcomponent = 20;       % specify the component(s) that should be plotted
% comp           = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'pca';
cfg.updatesens = 'no';
cfg.channel = 'meg';
comp = ft_componentanalysis(cfg, data);

% cfg           = [];
% cfg.component = 1:20;       % specify the component(s) that should be plotted
% cfg.layout    = lay;
% cfg.comment   = 'no';
% ft_topoplotIC(cfg, comp)

cfg = [];
cfg.viewmode = 'component';
cfg.layout = lay;
ft_databrowser(cfg, comp);

end