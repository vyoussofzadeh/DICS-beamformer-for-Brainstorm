% cd(dirnamesubj);

%% Do SPM inverse solutions
Ndip = 250;
allmeshvert_mni = D.inv{1}.mesh.tess_mni.vert;
Fs = D.fsample;
twinlength = 0.1; % sec
foi = [0 40]; % Hz
anatcol = [238 203 193]/255;

inv_conditions = {'Undefined'};
inv_modalities = {{'MEG'} {'MEG'} {'MEG'}};
inv_type =       {'MMN','COH','EBB'};

% EBB : Bayesian beamformer
inv_fboi = foi; % All since data already filtered

% tmax = max(D.time);
Sp = cell(length(inv_type),1);

for val = 1:length(inv_type)
    tonset = toi(2,1);
    D.inv{val}.inverse = [];   % Clear to be safe!
    D.val = val;
    D.inv{val}.comment = {sprintf('Ind: Val%d: Mod:%s Inv:%s',val,cat(2,inv_modalities{val}{:}),inv_type{val})};
    if (val > 1)
        D.inv{val}.mesh    = D.inv{1}.mesh;
        D.inv{val}.datareg = D.inv{1}.datareg;
        D.inv{val}.forward = D.inv{1}.forward;
        D.inv{val}.gainmat = D.inv{1}.gainmat;
    end
    D.inv{val}.inverse.trials   = inv_conditions;
    D.inv{val}.inverse.Han      = 1;
    %         D.inv{val}.inverse.Han      = 0;
    D.inv{val}.inverse.lpf      = inv_fboi(1);
    D.inv{val}.inverse.hpf      = inv_fboi(2);
    D.inv{val}.inverse.type     = inv_type{val};
    D.inv{val}.inverse.modality = inv_modalities{val};
    
    Ntr = 0;
    while(tonset + twinlength < toi(2,2))
        inv_twin = 1000*(tonset + [0 twinlength]);
        D.inv{val}.inverse.woi = inv_twin;
        D = spm_eeg_invert(D);
        tonset = tonset + twinlength;
        Jsol = D.inv{val}.inverse.J{1}*D.inv{val}.inverse.T';
        Ntr = Ntr + 1;
        if (Ntr == 1)
            Sp{val} = 0;
            NFFT = size(Jsol,2);
            freq = Fs/2*linspace(0,1,NFFT/2+1);
            % [~,indsel] = min(abs(freq - foi));
            indsel = find(freq >= foi(1) & freq <= foi(2));
        end
        disp(Ntr);
        tmp = fft(Jsol,NFFT,2);
        Sp{val} = Sp{val} + sum(abs(tmp(:,indsel)).^2,2);
    end
    Sp{val} = Sp{val}/(Ntr*length(indsel));
    %     figure(val); plot_bsurf_allviews(mesh.tess_ctx, anatcol, Sp{val});
    
end

%- SPM
figure
for val = 1:length(inv_type)
    Is = D.inv{val}.inverse.Is;
    Jabs = Sp{val};
    [~,ii] = sort(Jabs,'descend');
    ii = ii(1:Ndip);
    ivert_spm{val} = Is(ii);
    subplot(2,2,val); spm_mip(Jabs(ivert_spm{val}),allmeshvert_mni(ivert_spm{val},:)',6);
    axis image; colormap gray;
    title(sprintf('%s: avg \\omega in [%d,%d] Hz', inv_type{val}, foi));
end

%% LCMV-Beamformer
[individual_headmodel, sens] = ft_prepare_vol_sens(individual_headmodel, sens, 'channel', t_data.app.label);
method = 'lcmv';
cfg                  = [];
cfg.method           = 'lcmv';
cfg.grid             = individual_grid; % leadfield, which has the grid information
cfg.headmodel        = individual_headmodel; % volume conduction model (headmodel)
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.lambda      = '0%';
cfg.lcmv.keepfilter  = 'yes';
sourceAll = ft_sourceanalysis(cfg, t_data.app);
cfg.grid.filter = sourceAll.avg.filter;
s_data.bsl = ft_sourceanalysis(cfg, t_data.bsl);
s_data.pst = ft_sourceanalysis(cfg, t_data.pst);

s_data2.bsl  = ft_sourcedescriptives([], s_data.bsl); % to get the neural-activity-index
s_data2.pst  = ft_sourcedescriptives([], s_data.pst); % to get the neural-activity-index

cfg = [];
cfg.parameter = 'pow';
% cfg.operation = '(x1-x2)/(x1+x2)';
cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
source_ratio = ft_math(cfg,s_data.pst,s_data.bsl);
source_ratio.pow(source_ratio.pow>0)=0;
source_ratio.pow = abs(source_ratio.pow);

figure
m = source_ratio.pow;
bnd.pnt = individual_grid.pos;
bnd.tri = individual_grid.tri;
ft_plot_mesh(bnd, 'vertexcolor', m);
colorbar

clear savepath
savepath{1} = [method,'_left'];
savepath{2} = [method,'_right'];

%- prepare for plotting with spm_mip
% source_ratio.pow(source_ratio.pow>0)=0;
Jabs_lcmv = abs(source_ratio.pow);
view([180 0])
colormap('hot')
hcp_write_figure([savepath{1},'.png'], gcf, 'resolution', 300);
view([0 0])
hcp_write_figure([savepath{2},'.png'], gcf, 'resolution', 300);

[~,ii] = sort(Jabs_lcmv,'descend');
Ndip = 250;
ii = ii(1:Ndip);
Is = D.inv{1}.inverse.Is;
ivert_lcmv = Is(ii);

% x = coh;
% p = 70; %// percent of values that should become 1
% threshold = prctile(x,p);
% x_quant = x>=threshold;
% ivert_lcmv = x; 
% ivert_lcmv(~x_quant)=0;

% - Display results
figure, spm_mip(Jabs_lcmv(ivert_lcmv),allmeshvert_mni(ivert_lcmv,:)',6);
axis image; colormap gray;
% title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'LCMV', [1,30]));


%% Connectivity
addpath(genpath(ft_old));
addpath(genpath(fullfile(cd_org,'functions')));

%
mtd = 'plv';
mtd_par = 'plvspctrm';
% mtd = 'coh'; 'plvspctrm';

% mtd = 'coh';
% mtd_par = 'cohspctrm';

% mtd_par = 'cohspctrm';
gtm = 'eigenvector_cent';
% gtm = 'degrees';
% gtm = 'clustering_coef';

conn_par = [];
conn_par.method   = mtd;
conn_par.idx      = mtd_par;
conn_par.complex  = [];
conn_par.complex  = 'abs';
% conn_par.complex  = 'complex';

net_par.gtm       = gtm ; in2 = 2;
net_par.threshold = 0.7;
[conn_bsl, network_bsl] = vy_conn(s_data2.bsl,conn_par,net_par);
[conn_pst, network_pst] = vy_conn(s_data2.pst,conn_par,net_par);

ec_diff = [];
ec_diff.(gtm) = zeros(size(s_data.bsl.avg.pow,1),1);
ec_diff.pos     = s_data.bsl.pos;
ec_diff.inside  = s_data.bsl.inside;
ec_diff.(gtm)(ec_diff.inside==1) = zscore(network_pst.(gtm) - network_bsl.(gtm));
ec_diff.(gtm)(ec_diff.(gtm)<0)=0;

addpath(genpath(spm_path))

Jabs_net = abs(ec_diff.(gtm));

figure
m = ec_diff.(gtm);
bnd.pnt = individual_grid.pos;
bnd.tri = individual_grid.tri;
ft_plot_mesh(bnd, 'vertexcolor', m);
colorbar
title('conn analysis');

[~,ii] = sort(Jabs_net,'descend');
Ndip = 150;
ii = ii(1:Ndip);
Is = D.inv{1}.inverse.Is;
Jabs_net1 = Is(ii);
figure, spm_mip(Jabs_net(Jabs_net1),allmeshvert_mni(Jabs_net1,:)',6);
axis image; colormap gray;
title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'PLV-Network', [1,30]));


% rmpath(genpath(spm_path))
% cfg1 = [];
% % cfg1.outparam = 'powcorr'; cfg1.inparam = 'powcov'; cfg1.complex = 'absimag';
% cfg1.outparam = 'plvspctrm'; cfg1.inparam = 'crsspctrm'; cfg1.complex = 'abs';
% % cfg1.outparam = 'cohspctrm'; cfg1.inparam = 'crsspctrm'; cfg1.complex = 'abs';
% conn_bsl  = vy_connectivity(cfg1, s_data2.bsl);
% conn_pst  = vy_connectivity(cfg1, s_data2.pst);
% 
% % gtm = 'degrees';
% gtm = 'eigenvector_cent';
% % gtm = 'betweenness';
% % gtm = 'pagerank';
% % gtm = 'clustering_coef';
% % gtm = 'subgraph_centrality';
% % gtm = 'kcoreness_centrality_bu';
% cfg = [];
% cfg.method    = gtm;
% cfg.parameter = cfg1.outparam;
% cfg.threshold = 0.7;
% network_bsl = ft_networkanalysis(cfg,conn_bsl);
% network_pst = ft_networkanalysis(cfg,conn_pst);
% 
% cfg = [];
% cfg.parameter = gtm;
% % cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
% % cfg.operation = '(x1/x2)'; % sourceA divided by sourceB
% cfg.operation = '(x1-x2)'; % sourceA divided by sourceB
% % cfg.operation = 'x1-0*x2'; % sourceA divided by sourceB
% network_diff = ft_math(cfg, network_pst, network_bsl);
% network_diff.(gtm)(network_diff.(gtm)<0)=0;

%
% vs_bsl = cell2mat(s_data2.bsl.avg.mom);
% vs_pst = cell2mat(s_data2.pst.avg.mom);

% cor_bsl  = atan(corr((vs_bsl)'));
% cor_pst  = atan(corr((vs_pst)'));

% cor_bsl  = atan(abs(corr((vs_bsl)')));
% cor_pst  = atan(abs(corr((vs_pst)')));
% 
% adj_bsl = double(cor_bsl(:,:)>0.3);
% adj_pst = double(cor_pst(:,:)>0.3);
% 
% ec_bsl = eigenvector_centrality_und(adj_bsl);
% ec_pst = eigenvector_centrality_und(adj_pst);
% 
% gtm = 'eigenvector_cent';
% % gtm = 'degrees';
% 
% % figure, hist(ec_pst,10000)
% ec_diff = [];
% ec_diff.(gtm) = zeros(size(s_data.bsl.avg.pow,1),1);
% ec_diff.pos     = s_data.bsl.pos;
% ec_diff.inside  = s_data.bsl.inside;
% ec_diff.(gtm)(ec_diff.inside==1) = zscore(ec_pst - ec_bsl);
% ec_diff.(gtm)(ec_diff.(gtm)<0)=0;
% 
% addpath(genpath(spm_path))
% 
% Jabs_net = abs(ec_diff.(gtm));
% 
% figure
% m = ec_diff.(gtm);
% bnd.pnt = individual_grid.pos;
% bnd.tri = individual_grid.tri;
% ft_plot_mesh(bnd, 'vertexcolor', m);
% colorbar


% conn_par.conn_thre = 0.9;
% conn_par.idx = 'plvspctrm';
% conn_ratio = vy_connvis(conn_pst,conn_bsl,conn_par, individual_headmodel);
% view([156,47])

%%
% vs_bsl = cell2mat(s_data2.bsl.avg.mom);
% vs_pst = cell2mat(s_data2.pst.avg.mom);
% 
% figure,plot(vs')
% 
% % cor_bsl  = (corr(envelope(vs_bsl)'));
% % cor_pst  = (corr(envelope(vs_pst)'));
% 
% cor_bsl  = corr((vs_bsl)');
% cor_pst  = corr((vs_pst)');
% 
% cor_diff = cor_pst - cor_bsl; 
% figure, imagesc(cor_diff)
% ec = eigenvector_centrality_und(cor_diff);
% 
% ec_bsl = eigenvector_centrality_und(cor_bsl);
% ec_pst = eigenvector_centrality_und(cor_pst);
% 
% ec_diff = (ec_pst - ec_bsl);
% 
% Jabs_net = ec_diff;
% 
% [~,ii] = sort(Jabs_net,'descend');
% Ndip = 150;
% ii = ii(1:Ndip);
% Is = D.inv{1}.inverse.Is;
% Jabs_net1 = Is(ii);
% figure, spm_mip(Jabs_net(Jabs_net1),allmeshvert_mni(Jabs_net1,:)',6);
% axis image; colormap gray;

%% mne
cfg               = [];
cfg.method        = 'mne';
cfg.grid          = individual_grid;
cfg.vol           = individual_headmodel;
cfg.mne.prewhiten = 'yes';
cfg.mne.lambda    = 3;
cfg.mne.scalesourcecov = 'yes';
cfg.mne.keepfilter  = 'yes';
sourceAll = ft_sourceanalysis(cfg, t_data.app);
cfg.grid.filter = sourceAll.avg.filter;
s_data_mne.bsl = ft_sourceanalysis(cfg, t_data.bsl);
s_data_mne.pst = ft_sourceanalysis(cfg, t_data.pst);

% cfg = [];
% cfg.projectmom = 'yes';
% s_data_mne1.bsl = ft_sourcedescriptives(cfg,s_data_mne.bsl);
% s_data_mne1.pst = ft_sourcedescriptives(cfg, s_data_mne.pst);
% 
% sdDIFF = s_data_mne1.pst;
% sdDIFF.avg.pow = s_data_mne1.pst.avg.pow - s_data_mne1.bsl.avg.pow;
% sdDIFF.tri = individual_grid.tri;
% sdDIFF.pnt = individual_grid.pos;

% cfg = [];
% cfg.mask = 'avg.pow';
% % cfg.funparameter = 'avg.pow';
% ft_sourcemovie(cfg,sdDIFF);

s_data_mne.pst.avg.pow = mean(s_data_mne.pst.avg.pow,2);
s_data_mne.bsl.avg.pow = mean(s_data_mne.bsl.avg.pow,2);

cfg = [];
cfg.parameter = 'pow';
% cfg.operation = '(x1-x2)/(x2)';
cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
source_ratio = ft_math(cfg,s_data_mne.pst,s_data_mne.bsl);

figure
m = source_ratio.pow;
bnd.pnt = individual_grid.pos;
bnd.tri = individual_grid.tri;
ft_plot_mesh(bnd, 'vertexcolor', m);

%- prepare for plotting with spm_mip
source_ratio.pow(source_ratio.pow<0)=0;
Jabs_mne = (source_ratio.pow);

[~,ii] = sort(Jabs_mne,'descend');
Ndip = 250;
ii = ii(1:Ndip);
Is = D.inv{1}.inverse.Is;
ivert_mne = Is(ii);

% - Display results
figure,
in = ivert_mne;
spm_mip(Jabs_mne(in),allmeshvert_mni(in,:)',6);axis image; colormap gray;
axis image; colormap gray;
title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'mne', [0,30]));

%% DICS
freqs = {[6,13],[13 25],[6 30]};
for val = 1:length(freqs)
    
    cfg = [];
    cfg.headmodel  = individual_headmodel;  % from FT
    cfg.grad      = sens; % from FT
    cfg.senstype  = 'meg';
    cfg.grid      = individual_grid;
    cfg.method    = 'dics';
    cfg.dics.fixedori      = 'yes';
    %     cfg.dics.lambda = '30%';
    cfg.frequency = freqs{val};
    cfg.latency   = toi(2,:);
    sourceA = ft_sourceanalysis(cfg, w_data);
    cfg.latency   = [-0.380 -0.080];
    sourceB = ft_sourceanalysis(cfg, w_data);
    
    % FT_MATH requires the time axis to be the same
    sourceA.time = 0;
    sourceB.time = 0;
    
    cfg = [];
    cfg.parameter = 'pow';
    cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
    source_diff_dics = ft_math(cfg, sourceA, sourceB);
    
    Jabs_DICS{val} = abs(source_diff_dics.pow);
    
    [~,ii] = sort(Jabs_DICS{val},'descend');
    Ndip = 250;
    ii = ii(1:Ndip);
    Is = D.inv{1}.inverse.Is;
    ivert_dics{val} = Is(ii);
    
    %     figure(3), subplot(1,3,i); spm_mip(Jabs(ivert_dics{1}),allmeshvert_mni(ivert_dics{i},:)',6);
    %     title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'DICS', freqs{i}));
    %     axis image; colormap gray;
end

%% Display results, all
%- SPM
figure
clear ivert_spm
foi = [0 30]; % Hz
for val = 1:length(inv_type)
%     Is = D.inv{val}.inverse.Is;
    Jabs = Sp{val};
%     [~,ii] = sort(Jabs,'descend');
%     ii = ii(1:Ndip);
%     ivert_spm{val} = Is(ii);
    
    x = Jabs;
    p = 95; %// percent of values that should become 1
    threshold = prctile(x,p);
    x_quant = x>=threshold;
    x1 = x;
    x1(~x_quant)=0;
    subplot(2,2,val); spm_mip(x1(x_quant==1),allmeshvert_mni(x_quant==1,:)',6);
    
%     subplot(2,2,val); spm_mip(Jabs(ivert_spm{val}),allmeshvert_mni(ivert_spm{val},:)',6);
    axis image; colormap gray;
    title(sprintf('%s: avg \\omega in [%d,%d] Hz', inv_type{val}, foi));
end
savefig = fullfile(outputdir1,['spm_source_',subj,'_',run]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

%- LCMV
figure,
in = ivert_lcmv;
subplot(2,2,1), spm_mip(Jabs_lcmv(in),allmeshvert_mni(in,:)',6);axis image; colormap gray;
title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'LCMV', [1,40]));

% - DICS
for val = 1:length(freqs)
    Jabs = Jabs_DICS{val}; in = ivert_dics{val};
    subplot(2,2,val+1), spm_mip(Jabs(in),allmeshvert_mni(in,:)',6); axis image; colormap gray;
    title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'DICS', freqs{val}));
    axis image; colormap gray;
end
savefig = fullfile(outputdir1,['bf_',subj,'_',run]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

%- network
figure,
in = Jabs_net1;
spm_mip(Jabs_net(in),allmeshvert_mni(in,:)',6); axis image; colormap gray;
title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'network (PLV+EC)', [1,40]));
savefig = fullfile(outputdir1,['net_',subj,'_',run]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

% %- mne
figure,
in = ivert_mne;
spm_mip(Jabs_mne(in),allmeshvert_mni(in,:)',6);axis image; colormap gray;
title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'mne', [1,40]));
savefig = fullfile(outputdir1,['mne_',subj,'_',run]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

%% saving data
savepath = ['source_coh_',subj,'_',run,'.mat'];
coh = Sp{2};
save(savepath, 'coh', '-v7.3');

savepath = ['source_lcmv_',subj,'_',run,'.mat'];
lcmv = Jabs_lcmv;
save(savepath, 'lcmv', '-v7.3');

savepath = ['source_mne_',subj,'_',run,'.mat'];
mne = Jabs_mne;
save(savepath, 'mne', '-v7.3');
 
savepath = ['source_dics_',subj,'_',run,'.mat'];
dics = Jabs_DICS;
save(savepath, 'dics','freqs', '-v7.3');

savepath = ['source_network_',subj,'_',run,'.mat'];
network = Jabs_net;
save(savepath, 'network', '-v7.3');

%% Write result to SPM volumes
% D.inv{val}.contrast.smooth = 12;
% % D.inv{val}.contrast.rms = 1;
% % D.inv{val}.contrast.scalefactor = 1;
% D = spm_eeg_inv_Mesh2Voxels(D);
% SourceImgs{val}{ss} = strvcat(D.inv{val}.contrast.fname);

%% parcellation

%% PCC
% clear ivert_pcc
% freqs = {[6,13],[13 25],[6 30]};
% for val = 1:length(freqs)
%     
%     cfg = [];
%     cfg.headmodel  = individual_headmodel;  % from FT
%     cfg.grad      = sens; % from FT
%     cfg.senstype  = 'meg';
%     cfg.grid      = individual_grid;
%     cfg.method    = 'pcc';
%     cfg.pcc.fixedori      = 'yes';
%     cfg.keeptrials        = 'yes';
%     %     cfg.dics.lambda = '30%';
%     cfg.frequency = freqs{val};
%     cfg.latency   = toi(2,:);
%     sourceA = ft_sourceanalysis(cfg, w_data);
%     cfg.latency   = [-0.380 -0.080];
%     sourceB = ft_sourceanalysis(cfg, w_data);
%     
%     sourceA1 = ft_sourcedescriptives([], sourceA); % to get the neural-activity-index
%     sourceB1 = ft_sourcedescriptives([], sourceB); % to get the neural-activity-index
%     
%     %% compute connectivity
% %     cfg         = [];
% %     cfg.method  ='coh';
% %     cfg.complex = 'absimag';
% %     source_conn_A = ft_connectivityanalysis(cfg, sourceA);
% %     source_conn_B = ft_connectivityanalysis(cfg, sourceB);
%     
%     %%
%     
%     % FT_MATH requires the time axis to be the same
%     sourceA1.time = 0;
%     sourceB1.time = 0;
%     
%     cfg = [];
%     cfg.parameter = 'pow';
%     cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
%     source_diff_pcc = ft_math(cfg, sourceA1, sourceB1);
%     
%     Jabs_PCC{val} = abs(source_diff_pcc.pow);
%     
%     [~,ii] = sort(Jabs_DICS{val},'descend');
%     Ndip = 250;
%     ii = ii(1:Ndip);
%     Is = D.inv{1}.inverse.Is;
%     ivert_pcc{val} = Is(ii);
%     
%     %     figure(3), subplot(1,3,i); spm_mip(Jabs(ivert_dics{1}),allmeshvert_mni(ivert_dics{i},:)',6);
%     %     title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'DICS', freqs{i}));
%     %     axis image; colormap gray;
% end
% % - PCC
% figure,
% for val = 1:length(freqs)
%     Jabs = Jabs_PCC{val}; in = ivert_pcc{val};
%     subplot(2,2,val+1), spm_mip(Jabs(in),allmeshvert_mni(in,:)',6); axis image; colormap gray;
%     title(sprintf('%s: avg \\omega in [%d,%d] Hz', 'PCC', freqs{val}));
%     axis image; colormap gray;
% end
