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

%% reshaping vs
tvs.bsl = reshape(vs_tr.bsl.trial,size(vs_tr.bsl.trial,2), trl*size(vs_tr.bsl.trial,3));
tvs.pst = reshape(vs_tr.pst.trial,size(vs_tr.pst.trial,2), trl*size(vs_tr.pst.trial,3));

%%
cor_bsl  = (abs(corr((tvs.bsl)')));
cor_pst  = (abs(corr((tvs.pst)')));

ec_bsl = eigenvector_centrality_und(cor_bsl);
ec_pst = eigenvector_centrality_und(cor_pst);

% figure, hist(ec_pst,10000)

ec_diff = [];
ec_diff.eigenvector_cent = zeros(size(s_data.bsl.pos,1),1);
ec_diff.pos     = template_grid.pos;
ec_diff.dim     = template_grid.dim;
ec_diff.inside  = template_grid.inside;
ec_diff.eigenvector_cent(ec_diff.inside==1) = zscore(ec_pst - ec_bsl);
ec_diff.eigenvector_cent(ec_diff.eigenvector_cent<0)=0;

%%
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.taper      = 'hanning';
% cfg.taper      = 'dpss';
cfg.tapsmofrq  = 3;
% cfg.foilim     = [18 23];
cfg.foilim     = [1 40];
cfg.pad        = 3;
vs_tr_freq.bsl = ft_freqanalysis(cfg, vs_tr.bsl);
vs_tr_freq.pst = ft_freqanalysis(cfg, vs_tr.pst);

%%
mean_freq_vs.bsl = squeeze(mean(abs(squeeze(mean(vs_tr_freq.bsl.fourierspctrm,1))),1));
mean_freq_vs.pst = squeeze(mean(abs(squeeze(mean(vs_tr_freq.pst.fourierspctrm,1))),1));
figure,plot(vs_tr_freq.bsl.freq,mean_freq_vs.pst - mean_freq_vs.bsl);

%%
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.taper      = 'hanning';
% cfg.taper      = 'dpss';
cfg.tapsmofrq  = 3;
cfg.foilim     = [12 23];
% cfg.foilim     = [1 40];
cfg.pad        = 3;
vs_tr_freq1.bsl = ft_freqanalysis(cfg, vs_tr.bsl);
vs_tr_freq1.pst = ft_freqanalysis(cfg, vs_tr.pst);

%%
par = 'wplispctrm';
cfg           = [];
% % cfg.method    = 'coh';
% % cfg.complex   = 'absimag';
cfg.method        = 'wpli_debiased'; cfg.complex   = par;
cfvs_source_bsl   = ft_connectivityanalysis(cfg, vs_tr_freq1.bsl);
cfvs_source_act   = ft_connectivityanalysis(cfg, vs_tr_freq1.pst);

%%
par = 'wpli_debiasedspctrm';
% par = 'cohspctrm';

gtm = 'eigenvector_cent';
% gtm =  'degrees';

cfg           = [];
cfg.method    = gtm;
cfg.parameter = par;
cfg.threshold = .1;
network_bsl = ft_networkanalysis(cfg,cfvs_source_bsl);
network_pst = ft_networkanalysis(cfg,cfvs_source_act);

%%
cfg = [];
cfg.parameter = gtm;
cfg.operation = 'x1-x2';
network_diff = ft_math(cfg,network_pst,network_bsl);
network_diff.pos     = template_grid.pos;
network_diff.dim     = template_grid.dim;
network_diff.inside  = template_grid.inside;


%% visualize
network_diff1 = network_diff;
network_diff1.(gtm) = zeros(size(s_data.bsl.pos,1),1);
network_diff1.(gtm)(network_diff1.inside==1) = squeeze(mean(network_diff.(gtm),2));
network_diff1.(gtm)(network_diff1.(gtm)<0)=0;

param = [];
param.mask = gtm;
param.loc = 'max';
network_int = vy_source_plot(network_diff1,template_mri,param,2);

vy_mapvisualisation(network_int,gtm,0.3, []);

%%
% using older ver of ft for network analysis
restoredefaultpath
addpath(genpath(ft_old));
addpath(genpath('.\functions'));

%%
mtd = 'plv';
mtd_par = 'plvspctrm';

% mtd = 'coh';
% mtd_par = 'cohspctrm';
gtm = 'eigenvector_cent';

conn_par = [];
conn_par.method   = mtd;
conn_par.idx      = mtd_par;
conn_par.complex  = [];

net_par.gtm       = gtm ; in2 = 2;
net_par.threshold = 0.7;
[source_conn_bsl, network_bsl] = vy_conn(s_data2.bsl,conn_par,net_par);
[source_conn_pst, network_pst] = vy_conn(s_data2.pst,conn_par,net_par);
clear s_data2

%%
cfg = [];
cfg.parameter = gtm;
cfg.operation = 'x1-x2';
network_diff_lcmv = ft_math(cfg,network_pst,network_bsl);
network_diff_lcmv.pos     = template_grid.pos;
network_diff_lcmv.dim     = template_grid.dim;
network_diff_lcmv.inside  = template_grid.inside;
network_diff_lcmv.eigenvector_cent(network_diff_lcmv.eigenvector_cent<0)=0;

%% revert to the latest ft!
restoredefaultpath
addpath(genpath('.\functions'));
addpath .\Data_file;

ft_path = 'F:\My Matlab\My codes\My GitHub\fieldtrip_041718\fieldtrip-master';
addpath(ft_path);
ft_defaults % this loads the rest of the defaults;

hcp_path = 'F:\My Matlab\MEG\HCP\megconnectome-3.0';
addpath(genpath(hcp_path));

atlas_path = fullfile(ft_path,'template','atlas');
atlas = ft_read_atlas(fullfile(atlas_path,'aal\ROI_MNI_V4.nii'));

load temp_grid_8mm % from, vy_warping()
template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %

DestDirectory = 'H:\VNS'; % saving directory
%%
outputdir1 = fullfile(outputdir,'network');
if exist(outputdir1, 'file') == 0, mkdir(outputdir1), end
savedata = fullfile(outputdir1,['n_',subj,'_',run,'.mat']);
save(outputdir1, 'network_diff_lcmv', '-v7.3');

mtd = 'network_evc';

savepath = fullfile(outputdir1,[mtd,'_',subj,'_',run,'.mat']);
save(savepath, 'network_diff_lcmv', '-v7.3');

param = [];
param.mask = gtm;
param.loc = 'max';
network_int_lcmv = vy_source_plot(network_diff_lcmv,template_mri,param,2);
savefig = fullfile(outputdir1,[mtd,'_1_',subj,'_',run]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

clear savepath
savepath{1} = fullfile(outputdir1,[mtd,'_2_',subj,'_',run]);
savepath{2} = fullfile(outputdir1,[mtd,'_3_',subj,'_',run]);
vy_mapvisualisation(network_int_lcmv,gtm,0.3, savepath);

savenii = fullfile(outputdir1,['n_',subj,'_',run,'.nii']);
vy_savenifti(network_int_lcmv, gtm, savenii);

%% parcellation - aal (132 rois)
% network_diff_lcmv.eigenvector_cent = (network_diff_lcmv.eigenvector_cent)./max(network_diff_lcmv.eigenvector_cent);
% [~, data_intpar, coor] = vy_parcellate(network_diff_lcmv, atlas, gtm);
% data_intpar.eigenvector_centdimord = 'chan';
% 
% savepath = fullfile(outputdir1,[mtd,'_par',subj,'_',run,'.mat']);
% save(savepath, 'data_intpar', '-v7.3');
% 
% clear savepath
% savepath{1} = fullfile(outputdir1,[mtd,'_par_1',subj,'_',run]);
% savepath{2} = fullfile(outputdir1,[mtd,'_par_2',subj,'_',run]);
% vy_mapvisualisation(data_intpar,gtm,0.3, savepath);
% 
% % save nii
% savenii = fullfile(outputdir1,['n_par',subj,'_',run,'.nii']);
% vy_savenifti(data_intpar, gtm, savenii);

%% Conn
conn_par.conn_thre = 0.9;
conn_ratio = vy_connvis(source_conn_pst,source_conn_bsl,conn_par, individual_headmodel, network_diff_lcmv);
view([156,47])
savepath = fullfile(outputdir1,['conn_temp',subj,'_',run]);
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
savepath = fullfile(outputdir1,[mtd,'_conn_',subj,'_',run,'.mat']);
save(savepath, 'conn_ratio', 'individual_headmodel','-v7.3');

%% ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, gtm);
% disp(ROI_sel)
% savepath = fullfile(outputdir1,['n_ROIs_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
