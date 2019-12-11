function f_data = vy_preprocess5(datafile,task)

% cfg                         = [];
% cfg.dataset                 = datafile;
% % cfg.trialfun                = 'ft_trialfun_general'; % this is the default
% switch task
%     case 1
%         cfg.trialfun = 'trialfun_CRM'; % this is the default
%     case 2
%         cfg.trialfun = 'trialfun_VerbGenAud'; % this is the default
%     case 3
%         cfg.trialfun = 'trialfun_VergGenVis'; % this is the default
% end
% % cfg.trialdef.eventtype      = epoch_type;
% % cfg.trialdef.eventvalue     = Evnt_IDs; % the value of the stimulus trigger for fully incongruent (FIC).
% cfg.trialdef.prestim        = 1; % in seconds
% cfg.trialdef.poststim       = 2; % in seconds
% cfg.trialdef.prestimTime = cfg.trialdef.prestim;
% cfg.trialdef.poststimTime = cfg.trialdef.poststim;
% cfg = ft_definetrial(cfg);

cfg = [];
cfg.dataset = datafile;
switch task
    case 1
        cfg.trialfun = 'trialfun_CRM'; % this is the default
    case 2
        cfg.trialfun = 'trialfun_VerbGenAud'; % this is the default
    case 3
        cfg.trialfun = 'trialfun_VergGenVis'; % this is the default
end
cfg.trialdef.prestimTime    = 1; % in seconds
cfg.trialdef.poststimTime   = 2; % in seconds
cfg.trialdef.plotresults    = 'yes';
cfg = ft_definetrial(cfg);

cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.dftfilter = 'yes';
cfg.hpfiltord = 3;
cfg.dftfreq = [60 120 180];
cfg.hpfreq = 0.1;
cfg.lpfreq = 40;

% cfg.hpfreq = 5;
% cfg.lpfreq = 30;
cfg.channel = {'MEG'};
cfg.demean = 'yes';
cfg.baselinewindow = [-0.45 0.0];
% cfg.coilaccuracy = 0;
f_data = ft_preprocessing(cfg);

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

