function stat = vy_source_stat_montcarlo(s_data)


%%
cfg = [];
cfg.parameter        = 'pow';
% cfg.dim              = s_data.pst.dim;
cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.statistic        = 'depsamplesT';
% cfg.statistic         = 'indepsamplesT'; 
% cfg.correctm         = 'cluster';
cfg.correctm         = 'fdr';
cfg.clusteralpha     = 0.001;
% cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;

% cfg = [];
% cfg.method            = 'montecarlo';           % use the Monte Carlo Method to calculate the significance probability
% cfg.statistic         = 'indepsamplesT';        % use the independent samples T-statistic as a measure to evaluate the effect at the sample level
% cfg.correctm          = 'cluster';
% cfg.clusteralpha      = 0.05;                   % alpha level of the sample-specific test statistic that will be used for thresholding
% cfg.clustertail       = 0;
% cfg.clusterstatistic  = 'maxsum';               % test statistic that will be evaluated under the permutation distribution.
% cfg.tail              = 0;                      % -1, 1 or 0 (default = 0); one-sided or two-sided test
% cfg.correcttail       = 'prob';                 % the two-sided test implies that we do non-parametric two tests
% cfg.alpha             = 0.05;                   % alpha level of the permutation test
% cfg.numrandomization  = 1000;                   % number of draws from the permutation distribution
% cfg.design            = TFRwavelet.trialinfo(:,1)'; % design matrix, note the transpose
% cfg.ivar              = 1;                      % the index of the independent variable in the design matrix
% cfg.channel           = {'MEG1243'};
% cfg.neighbours        = [];                     % there are no spatial neighbours, only in time and frequency
 
ntrials                       = numel(s_data.bsl.trial);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = 1:ntrials;
design(2,ntrials+1:2*ntrials) = 1:ntrials;
% cfg_neighb.method    = 'distance';
% cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, timeAll);
 
cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
stat         = ft_sourcestatistics(cfg,s_data.pst,s_data.bsl);
% stat.pos     = individual_grid.pos;% keep positions for plotting later
% mri_aligned  = output.mri.mri_aligned;