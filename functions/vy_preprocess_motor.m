function [f_data, ecg_data] = vy_preprocess_motor(cfg_main)

cfg                         = [];
cfg.dataset                 = cfg_main.datafile;
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = cfg_main.epochtype;
cfg.trialdef.eventvalue     = cfg_main.eventid; % the value of the stimulus trigger for fully incongruent (FIC).
cfg.trialdef.prestim        = 1; % in seconds
cfg.trialdef.poststim       = 1; % in seconds
cfg = ft_definetrial(cfg);

cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.dftfilter = 'yes';
cfg.hpfiltord = 3;
% cfg.dftfreq = [60 120 180];
cfg.hpfreq = cfg_main.hpfreq;
cfg.lpfreq = cfg_main.lpfreq;
% cfg.hpfreq = 2;
% cfg.lpfreq = 30;
cfg.channel = {'MEG'};
% cfg.channel = {'MEGGRAD'};
cfg.demean = 'yes';
% cfg.baselinewindow = [-0.45 0.0];
f_data = ft_preprocessing(cfg);


% cfg.channel = {'ECG064'};
% ecg_data = ft_preprocessing(cfg);


%% notch filter
%     switch task
%         case 1 %'DefNam'
%             cfg        = [];
%             cfg.bsfilter  = 'yes';
%             cfg.bsfreq = [30,33];
%             f_data1 = ft_preprocessing(cfg,f_data);
%     end

%%
%   cfg.preproc.bpfilter    = 'yes'
%   cfg.preproc.bpfreq      = [110 140]
%   cfg.preproc.bpfiltord   =  8
%   cfg.preproc.bpfilttype  = 'but'
%   cfg.preproc.rectify     = 'yes'
%   cfg.preproc.boxcar      = 0.2

% [cfg, artifact] = ft_artifact_jump(cfg);        % id jumps
% f_data = ft_rejectartifact(cfg, f_data);        % remove scanner (MEG) jumps

% cfg = [];
% cfg.dataset  = datafile;
% cfg.hpfilter = 'yes';
% cfg.lpfilter = 'yes';
% cfg.dftfilter = 'yes';
% cfg.hpfiltord = 3;
% cfg.dftfreq = [60 120 180];
% cfg.hpfreq = .1;
% cfg.lpfreq = 40;
% cfg.channel = {'MEG'};
% cfg.demean = 'yes';
% % cfg.baselinewindow = [-0.2 0.0];
% data = ft_preprocessing(cfg);
%
% cfg = [];
% cfg.dataset                 = datafile;
% cfg.trialfun                = 'ft_trialfun_general'; % this is the default
% cfg.trialdef.eventtype      = epoch_type;
% cfg.trialdef.eventvalue     = Evnt_IDs; % the value of the stimulus trigger for fully incongruent (FIC).
% cfg.trialdef.prestim        = 1; % in seconds
% cfg.trialdef.poststim       = 2; % in seconds
% cfg = ft_definetrial(cfg);
%
% trl = cfg.trl;
% cfg = [];
% cfg.trl = trl;
% f_data = ft_redefinetrial(cfg,data);
%
%
% cfg = [];
% cfg.dataset = datafile;
% cfg.trl = trl;
% [cfg, artifact] = ft_artifact_jump(cfg);        % id jumps
% f_data = ft_rejectartifact(cfg, f_data);        % remove scanner (MEG) jumps

