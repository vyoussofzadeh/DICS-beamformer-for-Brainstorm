
% check if lcmv has done!
if ~exist('s_data2','var')
    [s_data, s_data2] = vy_source(t_data, individual_grid, individual_headmodel);
end

%% virtual sens (pre, post)
clear vs
s_bsl = ft_checkdata(s_data2.bsl, 'datatype', {'freqmvar' 'freq' 'source'});
v = cell2mat(s_bsl.mom);
trl = numel(s_data2.bsl);
vs.bsl.trial = reshape(v,trl,size(v,1)/trl,size(v,2));
vs.bsl.time = s_data2.bsl.time;
% vs.bsl.label = cellstr(num2str(individual_grid.pos(individual_grid.inside,:)));

s_pst = ft_checkdata(s_data2.pst, 'datatype', {'freqmvar' 'freq' 'source'});
v = cell2mat(s_pst.mom);
trl = numel(s_data2.pst);
vs.pst.trial = reshape(v,trl,size(v,1)/trl,size(v,2));
vs.pst.time = s_data2.pst.time;
% vs.pst.label = cellstr(num2str(individual_grid.pos(individual_grid.inside,:)));

%%
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.taper      = 'hanning';
% cfg.taper      = 'dpss';
cfg.tapsmofrq  = 2;
cfg.foilim     = [18 23];
% cfg.foilim     = [1 40];
cfg.pad        = 3;
vs_f.bsl = ft_freqanalysis(cfg, vs.bsl);
vs_f.pst = ft_freqanalysis(cfg, vs.pst);

%%
cfg           = [];
cfg.method    = 'coh'; 
% cfg.method    = 'plv'; 
% cfg.method        = 'wpli_debiased'; 

cfg.complex   = 'absimag';
vs_f_conn.bsl   = ft_connectivityanalysis(cfg, vs_f.bsl);
vs_f_conn.pst   = ft_connectivityanalysis(cfg, vs_f.pst);

%%
vs_f_conn.bsl.conn = squeeze(mean(vs_f_conn.bsl.cohspctrm,3));
vs_f_conn.pst.conn = squeeze(mean(vs_f_conn.pst.cohspctrm,3));

% vs_f_conn.bsl.conn = squeeze(mean(vs_f_conn.bsl.plvspctrm,3));
% vs_f_conn.pst.conn = squeeze(mean(vs_f_conn.pst.plvspctrm,3));

% vs_f_conn.bsl.conn = squeeze(mean(vs_f_conn.bsl.wpli_debiasedspctrm,3));
% vs_f_conn.pst.conn = squeeze(mean(vs_f_conn.pst.wpli_debiasedspctrm,3));

%% Conn visualization
figure,
% imagesc(vs_f_conn.pst.cohspctrm - vs_f_conn.bsl.cohspctrm);
imagesc(vs_f_conn.pst.conn - vs_f_conn.bsl.conn);
title('Averaged Conn')

%%
% aedge = conn_ratio.cohspctrm;
aedge = vs_f_conn.pst.conn - vs_f_conn.bsl.conn;

aedge(isnan(aedge))=0;
tedge = (aedge.* double(aedge > 0.96.*max(aedge(:))));

figure,
ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none'); alpha 0.5; camlight;
view ([-196 56])

ROI = individual_grid.pos(individual_grid.inside,:);
hold on
for i = 1:length(tedge)
    for j = 1:length(tedge)
        if tedge(i,j)> 0
            p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
            p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
            pts = [p1; p2];
            line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
        end
    end
end
view([90 0])
% title(tit)

%%
% using older ver of ft for network analysis
restoredefaultpath
addpath(genpath(ft_old));
addpath(genpath('.\functions'));

%%
mtd = 'plv';
idx = 'conn';
gtm = 'eigenvector_cent';
net_par.threshold = 0.7;

cfg = [];
cfg.method    = gtm;
cfg.parameter = idx;
cfg.threshold = net_par.threshold;
net_bsl = ft_networkanalysis(cfg,vs_f_conn.bsl);
net_pst = ft_networkanalysis(cfg,vs_f_conn.pst);

%%
cfg = [];
cfg.parameter = gtm;
cfg.operation = 'x1-x2';
network_diff_lcmv = ft_math(cfg,net_pst,net_bsl);
network_diff_lcmv.pos     = template_grid.pos;
network_diff_lcmv.dim     = template_grid.dim;
network_diff_lcmv.inside  = template_grid.inside;
network_diff_lcmv.eigenvector_cent(network_diff_lcmv.eigenvector_cent<0)=0;
network_diff_lcmv.dimord = 'pos';

a = find(template_grid.inside==1);
evc = zeros(length(template_grid.pos),1);
evc(a)= network_diff_lcmv.eigenvector_cent;
network_diff_lcmv.eigenvector_cent = evc;

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
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = param.mask;
cfg.location     = param.loc;
ft_sourceplot(cfg,network_diff_lcmv);
cfg.funparameter = 'pow';
ft_sourceplot(cfg,s_data.bsl);

param = [];
param.mask = gtm;
param.loc = 'max';
network_int_lcmv = vy_source_plot(network_diff_lcmv,template_mri,param,2);

vy_mapvisualisation(network_int_lcmv,gtm,0.3, []);
