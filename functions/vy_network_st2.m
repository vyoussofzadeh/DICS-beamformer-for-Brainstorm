BCT_path = 'F:\My Matlab\Network analysis\BTC_Brain Connectivity Toolbox\BCT\2016_01_16_BCT';
addpath(genpath(BCT_path))

%% source- single trial
s_data = vy_source_stat(t_data, individual_grid, individual_headmodel);

%%
% check if lcmv has done!
% if ~exist('s_data2','var')
%     [s_data, s_data2] = vy_source(t_data, individual_grid, individual_headmodel);
% end

%% virtual sens (pre, post)
vs_tr = [];
source1 = ft_checkdata(s_data.bsl, 'datatype', {'freqmvar' 'freq' 'source'});
vs = cell2mat(source1.mom);
trl = numel(s_data.bsl.trial);
vs_tr.bsl.trial = reshape(vs,trl,size(vs,1)/trl,size(vs,3));
vs_tr.bsl.time = s_data.bsl.time;
vs_tr.bsl.label = cellstr(num2str(individual_grid.pos(individual_grid.inside,:)));

source1 = ft_checkdata(s_data.pst, 'datatype', {'freqmvar' 'freq' 'source'});
vs = cell2mat(source1.mom);
trl = numel(s_data.bsl.trial);
vs_tr.pst.trial = reshape(vs,trl,size(vs,1)/trl,size(vs,3));
vs_tr.pst.time = s_data.bsl.time;
vs_tr.pst.label = cellstr(num2str(individual_grid.pos(individual_grid.inside,:)));

% vs_source_act = [];
% for i = 1:length(source_act.trial) % i = number trials
%     for x = 1:length(source_act.trial(i).mom) % x = number nodes
%         node(x,:) = source_act.trial(i).mom{x};
%     end
%     vs_source_act.trial{i} = node;
%     vs_source_act.time{1,i} = source_act.time;
% end
% vs_source_act.label = cellstr(num2str(headmodel_time.pos(headmodel_time.inside,:)));

%%
clear EC_pre EC_pst
for i=1:size(vs_tr.bsl.trial,1)
    disp([num2str(i),'/',num2str(size(vs_tr.bsl.trial,1))])
%     EC_pre(i,:) = eigenvector_centrality_und(plv_measure(squeeze(vs_tr.bsl.trial(i,:,:))));
%     EC_pst(i,:) = eigenvector_centrality_und(plv_measure(squeeze(vs_tr.pst.trial(i,:,:))));
    
    EC_pre(i,:) = eigenvector_centrality_und(corr(squeeze(vs_tr.bsl.trial(i,:,:))'));
    EC_pst(i,:) = eigenvector_centrality_und(corr(squeeze(vs_tr.pst.trial(i,:,:))'));
end

%%
sourcepreS1_d = s_data.bsl;
sourcepstS1_d = s_data.pst;
for i = 1:length(s_data.bsl.trial) % i = number trials
    sourcepreS1_d.trial(i).pow(sourcepreS1_d.inside) = (EC_pre(i,:))';
    sourcepstS1_d.trial(i).pow(sourcepstS1_d.inside) = (EC_pst(i,:))';
end

%% FDR correction
cfg = [];
cfg.parameter        = 'pow';
cfg.dim              = individual_grid.dim;
cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.statistic        = 'depsamplesT';
% cfg.statistic         = 'indepsamplesT'; 
cfg.correctm         = 'max';
% cfg.correctm         = 'fdr';
cfg.clusteralpha     = 0.001;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;
%
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
stat         = ft_sourcestatistics(cfg,sourcepreS1_d,sourcepstS1_d);
stat.pos     = template_grid.pos;% keep positions for plotting later
% mri_aligned  = output.mri.mri_aligned;
stat.inside = template_grid.inside;

%%
stat.stat(stat.stat < 0 )=0;

param = [];
param.mask = 'stat';
param.loc = 'max';
network_stat = vy_source_plot(stat,template_mri,param,2);

vy_mapvisualisation(network_stat,'stat',0.4, [],0);
