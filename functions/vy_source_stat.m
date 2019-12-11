function s_data = vy_source_stat(data, grid, vol)

% create spatial filter using the lcmv beamformer
% cfg                  = [];
% cfg.method           = 'lcmv';
% cfg.grid             = grid; % leadfield, which has the grid information
% cfg.headmodel         = vol; % volume conduction model (headmodel)
% cfg.keepfilter       = 'yes';
% cfg.lcmv.keepfilter  = 'yes';
% cfg.keeptrials       = 'yes';
% cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
% cfg.lcmv.lambda      = '5%';
% % cfg.senstype = 'meg';
% 
% sourceAll = ft_sourceanalysis(cfg, data.app);
% 
% cfg.grid.filter = sourceAll.avg.filter;
% cfg.rawtrial    = 'yes';
% s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
% s_data.pst      = ft_sourceanalysis(cfg, data.pst);
% s_data2.bsl     = ft_sourcedescriptives([], s_data.bsl); % to get the neural-activity-index
% s_data2.pst     = ft_sourcedescriptives([], s_data.pst); % to get the neural-activity-index


cfg = [];
cfg.method = 'lcmv';
cfg.lcmv.lambda = '5%';
cfg.grid = grid;
cfg.vol = vol;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.channel = data.bsl.label;
sourceavg = ft_sourceanalysis(cfg, data.app);

cfg = [];
cfg.method = 'lcmv';
cfg.lcmv.lambda = '0.1%';
cfg.grid = grid;
cfg.grid.filter = sourceavg.avg.filter;
cfg.rawtrial = 'yes';
cfg.vol = vol;
s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
s_data.pst      = ft_sourceanalysis(cfg, data.pst);
