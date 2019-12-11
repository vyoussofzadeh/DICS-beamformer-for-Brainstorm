
%% Borowsing data
cfg = [];
cfg.viewmode = 'vertical';
cfg.continuous = 'no';
ft_databrowser(cfg,raw_data);


%% Preprocessing data
%
cln_data_BAK = cln_data;
% cln_data = cln_data_BAK;
cfg              = [];
cfg.hpfilter     = 'yes';        % enable high-pass filtering
cfg.lpfilter     = 'yes';        % enable low-pass filtering
cfg.hpfreq       = 1;           % set up the frequency for high-pass filter
cfg.lpfreq       = 8;          % set up the frequency for low-pass filter
fcln_data        = ft_preprocessing(cfg,cln_data_BAK);
