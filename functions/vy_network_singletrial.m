BCT_path = 'F:\My Matlab\Network analysis\BTC_Brain Connectivity Toolbox\BCT\2016_01_16_BCT';
addpath(genpath(BCT_path))


%%
[s_data2] = vy_source_stat(t_data, individual_grid, individual_headmodel);

%% virtual sens (pre, post)
in = s_data2.bsl;
vs_tr = [];
source = ft_checkdata(in, 'datatype', {'freqmvar' 'freq' 'source'});
vs = cell2mat(source.mom);
trl = numel(in.trial);
vs1.trial = reshape(vs,trl,size(vs,1)/trl,size(vs,3));
vs1.time = in.time;
vs1.label = cellstr(num2str(individual_grid.pos(individual_grid.inside,:)));
tvs.bsl = vs1;

in = s_data2.pst;
vs_tr = [];
source = ft_checkdata(in, 'datatype', {'freqmvar' 'freq' 'source'});
vs = cell2mat(source.mom);
trl = numel(in.trial);
vs1.trial = reshape(vs,trl,size(vs,1)/trl,size(vs,3));
vs1.time = in.time;
vs1.label = cellstr(num2str(individual_grid.pos(individual_grid.inside,:)));
tvs.pst = vs1;

%%
cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
cfg.vartrllength     = 0;
cfg.covariance       = 'yes';
cfg.keeptrials       = 'yes';
time_vs_bsl         = ft_timelockanalysis(cfg, tvs.bsl);
time_vs_act         = ft_timelockanalysis(cfg, tvs.pst);

%% GC
% GCCA_path = 'F:\My Matlab\Granger Causality\GCCA3';
% addpath(genpath(fullfile(GCCA_path)));
% 
% [bic,aic] = cca_find_model_order(time_vs_bsl.trial,2,15);
% 
% Nr = size(time_vs_bsl.trial,1);
% Nl = size(time_vs_bsl.trial,2);
% [bic,aic] = cca_find_model_order_mtrial(time_vs_bsl.trial,Nr,Nl,2,15);
% disp(['best model order by Bayesian Information Criterion = ',num2str(bic)]);
% disp(['best model order by Aikaike Information Criterion = ',num2str(aic)]);
% nlags = max(aic,bic);
% 
% m = mod(Nl,nBwin);
% X = X(:,1:end-Nr*m);
% N = length(X);  % length of whole eeg data
% Nl = N/Nr;  % length of each realization
% 
% nBoot = 3;
% nBwin = 20;
% pval = 0.05;
% CORRTYPE = 1;   % Bonferroni
% ret = cca_partialgc_doi_bstrap(time_vs_bsl.trial,Nr,Nl,nlags,nBoot,nBwin,pval,CORRTYPE,1);

%% MVAR analysis
% mvar
% cfg             = [];
% cfg.order       = 5;
% cfg.toolbox     = 'bsmart';
% mdata           = ft_mvaranalysis(cfg, time_vs_bsl);

%% Freq analysis
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
% cfg.taper      = 'hanning';
cfg.taper      = 'dpss';
cfg.tapsmofrq  = 4;
cfg.foilim     = [1 20];
cfg.pad        = 2;
fvs_source_bsl = ft_freqanalysis(cfg, time_vs_bsl);
fvs_source_act = ft_freqanalysis(cfg, time_vs_act);

%% Conn analysis
cfg           = [];
% cfg.method    = 'coh'; cfg.complex   = 'imag'; mm = 'cohspctrm';
cfg.method        = 'wpli_debiased'; mm = 'wpli_debiasedspctrm';
cfvs_source_bsl   = ft_connectivityanalysis(cfg, fvs_source_bsl);
cfvs_source_act   = ft_connectivityanalysis(cfg, fvs_source_act);

%%
% cfvs_source_bsl.(mm) = abs(cfvs_source_bsl.(mm));
% cfvs_source_act.(mm) = abs(cfvs_source_act.(mm));

mcfvs_source_bsl = cfvs_source_bsl;
mcfvs_source_bsl.(mm) = squeeze(mean(mcfvs_source_bsl.(mm),3));

mcfvs_source_act = cfvs_source_act;
mcfvs_source_act.(mm) = squeeze(mean(mcfvs_source_act.(mm),3));

figure,
imagesc(mcfvs_source_act.(mm) - mcfvs_source_bsl.(mm));colorbar

%%
mcfvs_source_bsl.(mm)(isnan(mcfvs_source_bsl.(mm))) = 0;
mcfvs_source_act.(mm)(isnan(mcfvs_source_act.(mm))) = 0;

mcfvs_source_bsl.(mm) = abs(mcfvs_source_bsl.(mm));
mcfvs_source_act.(mm) = abs(mcfvs_source_act.(mm));

%%
% cfg            = [];
% cfg.operation  = '(x1-x2)';
% cfg.parameter  = 'cohspctrm';
% conn_ratio     = ft_math(cfg, mcfvs_source_act, mcfvs_source_bsl);

%%
ec_bsl = eigenvector_centrality_und(mcfvs_source_bsl.(mm));
ec_pst = eigenvector_centrality_und(mcfvs_source_act.(mm));

% figure, hist(ec_pst,10000)
rr  =s_data2.bsl.trial.pow;
ec_diff = [];
ec_diff.eigenvector_cent = zeros(size(rr,1),1);
ec_diff.pos     = template_grid.pos;
ec_diff.dim     = template_grid.dim;
ec_diff.inside  = template_grid.inside;
ec_diff.eigenvector_cent(ec_diff.inside==1) = zscore(ec_pst - ec_bsl);
ec_diff.eigenvector_cent(ec_diff.eigenvector_cent<0)=0;

gtm = 'eigenvector_cent';
param = [];
param.mask = gtm;
param.loc = 'max';
network_int_lcmv = vy_source_plot(ec_diff,template_mri,param,2);

vy_mapvisualisation(network_int_lcmv,gtm,0.5, []);
