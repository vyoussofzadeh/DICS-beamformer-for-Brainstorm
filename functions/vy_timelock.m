function t_data = vy_timelock(data)

%
cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
cfg.preproc.demean   = 'yes';    % enable demean to remove mean value from each single trial
cfg.keeptrials       = 'yes';
% cfg.vartrllength = 0;

t_data.pst            = ft_timelockanalysis(cfg, data.pst);
t_data.bsl            = ft_timelockanalysis(cfg, data.bsl);
t_data.all            = ft_timelockanalysis(cfg, data.all);

cfg.vartrllength     = 2;
t_data.app           = ft_timelockanalysis(cfg, data.app);
