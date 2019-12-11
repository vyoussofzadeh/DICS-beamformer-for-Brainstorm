function a_data = vy_ave(data)

tt = data.time{1};
idx = tt==0;

cfg                   = [];
cfg.covariance        = 'yes';
cfg.covariancewindow  = 'all';
cfg.preproc.baselinewindow = [tt(1),tt(idx)];
cfg.preproc.demean    = 'yes';    % enable demean to remove mean value from each single trial

% a_data.pst         = ft_timelockanalysis(cfg, data.pst);
% a_data.pst.elec    = elec;

% a_data.bsl         = ft_timelockanalysis(cfg, data.bsl);
% a_data.bsl.elec    = elec;

a_data         = ft_timelockanalysis(cfg, data);
% a_data.all.elec    = elec;

% a_data.aud         = ft_timelockanalysis(cfg, data.aud);
% a_data.aud.elec    = elec;

%%
% cfg               = [];
% % cfg.reref         = 'yes';
% % cfg.refchannel    = 'all';
% cfg.refmethod     = 'avg';
% a_data.all        = ft_preprocessing(cfg,a_data.all); 
% a_data.pst        = ft_preprocessing(cfg,a_data.pst); 
% a_data.bsl        = ft_preprocessing(cfg,a_data.bsl); 
% a_data.aud        = ft_preprocessing(cfg,a_data.aud); 
