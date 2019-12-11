BCT_path = 'F:\My Matlab\Network analysis\BTC_Brain Connectivity Toolbox\BCT\2016_01_16_BCT';
addpath(genpath(BCT_path))

% check if lcmv exists!
% if ~exist('s_data2','var')
[s_data, s_data2] = vy_source(t_data, individual_grid, individual_headmodel);
% end
%%
vs_bsl = cell2mat(s_data2.bsl.avg.mom);
vs_pst = cell2mat(s_data2.pst.avg.mom);

%%
% [u,s,v] = svd(vs_pst,'econ');
% vs_pst1 = u(:,1)'*vs_pst;
% 
% figure,plot(u')


%%
% addpath('F:\My Matlab\MEG\MEG-ROI-nets-master\MEG-ROI-nets-master');
% slvs_bsl = ROInets.remove_source_leakage(vs_bsl, 'householder');



% hvs_bsl = hilbert(vs_bsl);
% hvs_pst = hilbert(vs_pst);
% 
% [u1,v1,d1] = svd(hvs_bsl,'econ');
% [u2,v2,d2] = svd(hvs_pst,'econ');
% 
% ahvs_bsl = abs(hvs_bsl);
% ahvs_pst = abs(hvs_pst);
% 
% figure,plot(ahvs_bsl')
% 
% iSeed = 1;
% HB1 = imag(bsxfun(@times, hvs_bsl, conj(u1(iSeed,:))./abs(u1(iSeed,:))));
% HBo1 = -imag(HB1 .* ((1i*u1)./abs(u1)));
% % figure,plot(HBo1')
% 
% HB2 = imag(bsxfun(@times, hvs_pst, conj(u2(iSeed,:))./abs(u2(iSeed,:))));
% HBo2 = -imag(HB2 .* ((1i*u2)./abs(u2)));
% % figure,plot(HBo2')

%%
% % ref: https://github.com/brainstorm-tools/brainstorm3/blob/4ccb8490749efdb98ae2c87084f49059e1a23987/toolbox/connectivity/bst_connectivity.m#L447
% HA = hvs_bsl;
% HB = hvs_pst;
% 
% for node = 1:size(HA,1)
%     % Orthogonalize complex coefficients, based on Hipp et al. 2012
%     % HBo is the amplitude of the component orthogonal to HA
%     HBo = imag(bsxfun(@times, HB, conj(HA(node,:))./abs(HA(node,:))));
%     % The orthogonalized signal can be computed like this (not necessary here):
%     % HBos = real(HBo .* ((1i*HA)./abs(HA)));
%     % avoid rounding errors
%     HBo(abs(HBo./abs(HB))<2*eps)=0;
%     % Compute correlation coefficients
%     R_bsl(node,:) = correlate_dims(abs(HBo), abs(HA(node,:)), 2);
% end
% 
% for node = 1:size(HB,1)
%     % Orthogonalize complex coefficients, based on Hipp et al. 2012
%     % HBo is the amplitude of the component orthogonal to HA
%     HBo = imag(bsxfun(@times, HB, conj(HB(node,:))./abs(HB(node,:))));
%     % The orthogonalized signal can be computed like this (not necessary here):
%     % HBos = real(HBo .* ((1i*HA)./abs(HA)));
%     % avoid rounding errors
%     HBo(abs(HBo./abs(HB))<2*eps)=0;
%     % Compute correlation coefficients
%     R_pst(node,:) = correlate_dims(abs(HBo), abs(HB(node,:)), 2);
% end

%%
plv_bsl  = atan(plv_measure((vs_bsl)));
plv_pst  = atan(plv_measure((vs_pst)));

ec_bsl = eigenvector_centrality_und(plv_bsl);
ec_pst = eigenvector_centrality_und(plv_pst);
ec_d = eigenvector_centrality_und(plv_pst - plv_bsl);

%%
opt = [];
opt.method = 'partialcorr';
opt.type = 'spearman';
cor_bsl  = atan(correlate((vs_bsl)',opt));
cor_pst  = atan(correlate((vs_pst)',opt));

ec_bsl = eigenvector_centrality_und(cor_bsl);
ec_pst = eigenvector_centrality_und(cor_pst);


% figure,imagesc(plv_pst)

% cor_bsl  = atan(corr((vs_bsl)'));
% cor_pst  = atan(corr((vs_pst)'));

% cor_bsl  = atan(corr((R_bsl)'));
% cor_pst  = atan(corr((R_pst)'));

% cor_bsl  = atan(corr((vs_bsl)'));
% cor_pst  = atan(corr((vs_pst)'));
% 
% cor_bsl  = atan(corr((ahvs_bsl)'));
% cor_pst  = atan(corr((ahvs_pst)'));

% cor_bsl  = ((corr((vs_bsl)')));
% cor_pst  = ((corr((vs_pst)')));

% ec_bsl = eigenvector_centrality_und(cor_bsl);
% ec_pst = eigenvector_centrality_und(cor_pst);

% figure, hist(ec_pst,10000)

ec_diff = [];
ec_diff.eigenvector_cent = zeros(size(s_data.bsl.avg.pow,1),1);
ec_diff.pos     = template_grid.pos;
ec_diff.dim     = template_grid.dim;
ec_diff.inside  = template_grid.inside;
% ec_diff.eigenvector_cent(ec_diff.inside==1) = ec_d;
ec_diff.eigenvector_cent(ec_diff.inside==1) = zscore(ec_pst - ec_bsl);
ec_diff.eigenvector_cent(ec_diff.eigenvector_cent<0)=0;

%%

% % check if lcmv exists!
% % if ~exist('s_data2','var')
% % [s_data, s_data2] = vy_source(t_data, individual_grid, individual_headmodel);
% [s_data2] = vy_source_stat(t_data, individual_grid, individual_headmodel);
% % end
% 
% %% virtual sens (pre, post)
% in = s_data2.bsl;
% vs_tr = [];
% source = ft_checkdata(in, 'datatype', {'freqmvar' 'freq' 'source'});
% vs = cell2mat(source.mom);
% trl = numel(in.trial);
% vs1.trial = reshape(vs,trl,size(vs,1)/trl,size(vs,3));
% vs1.time = in.time;
% vs1.label = cellstr(num2str(individual_grid.pos(individual_grid.inside,:)));
% tvs.bsl = vs1;
% 
% in = s_data2.pst;
% vs_tr = [];
% source = ft_checkdata(in, 'datatype', {'freqmvar' 'freq' 'source'});
% vs = cell2mat(source.mom);
% trl = numel(in.trial);
% vs1.trial = reshape(vs,trl,size(vs,1)/trl,size(vs,3));
% vs1.time = in.time;
% vs1.label = cellstr(num2str(individual_grid.pos(individual_grid.inside,:)));
% tvs.pst = vs1;
% 
% %%
% tvs.bsl.trial = reshape(tvs.bsl.trial,size(tvs.bsl.trial,2), trl*size(tvs.bsl.trial,3));
% tvs.pst.trial = reshape(tvs.pst.trial,size(tvs.pst.trial,2), trl*size(tvs.pst.trial,3));
% 
% %%
% % vs_bsl = cell2mat(s_data2.bsl.avg.mom);
% % vs_pst = cell2mat(s_data2.pst.avg.mom);
% 
% %%
% % addpath('F:\My Matlab\MEG\MEG-ROI-nets-master\MEG-ROI-nets-master');
% % slvs_bsl = ROInets.remove_source_leakage(vs_bsl, 'householder');
% hvs_bsl = hilbert(tvs.bsl.trial);
% hvs_pst = hilbert(tvs.pst.trial);
% 
% ahvs_bsl = abs(hvs_bsl);
% ahvs_pst = abs(hvs_pst);
% 
% figure,plot(ahvs_bsl')
% 
% HB1 = imag(bsxfun(@times, hvs_bsl, conj(hvs_bsl)./abs(hvs_bsl)));
% % HB1(abs(HB1./abs(hvs_bsl))<2*eps)=0;
% HBo1 = -imag(HB1 .* ((1i*hvs_bsl)./abs(hvs_bsl)));
% % figure,plot(HBo1')
% 
% HB2 = imag(bsxfun(@times, hvs_pst, conj(hvs_pst)./abs(hvs_pst)));
% % HB2(abs(HB2./abs(hvs_pst))<2*eps)=0;
% HBo2 = -imag(HB2 .* ((1i*hvs_pst)./abs(hvs_pst)));
% % figure,plot(HBo2')
% 
% %%
% cor_bsl  = atan(corr((HB1)'));
% cor_pst  = atan(corr((HB2)'));
% 
% % cor_bsl  = atan(abs(corr((HBo1)')));
% % cor_pst  = atan(abs(corr((HBo2)')));
% 
% % cor_bsl  = atan(corr((ahvs_bsl)'));
% % cor_pst  = atan(corr((ahvs_pst)'));
% 
% cor_bsl  = atan(abs(corr((tvs.bsl.trial)')));
% cor_pst  = atan(abs(corr((tvs.pst.trial)')));
% 
% ec_bsl = eigenvector_centrality_und(cor_bsl);
% ec_pst = eigenvector_centrality_und(cor_pst);
% 
% % figure, hist(ec_pst,10000)
% 
% ec_diff = [];
% ec_diff.eigenvector_cent = zeros(size(s_data.bsl.avg.pow,1),1);
% ec_diff.pos     = template_grid.pos;
% ec_diff.dim     = template_grid.dim;
% ec_diff.inside  = template_grid.inside;
% ec_diff.eigenvector_cent(ec_diff.inside==1) = (ec_pst - ec_bsl);
% ec_diff.eigenvector_cent(ec_diff.eigenvector_cent<0)=0;

%%
outputdir1 = fullfile(outputdir,'network');
if exist(outputdir1, 'file') == 0, mkdir(outputdir1), end
savedata = fullfile(outputdir1,['n_',subj,'_',run,'.mat']);
% save(outputdir1, 'ec_diff', '-v7.3');

mtd = 'network_evc';

savepath = fullfile(outputdir1,[mtd,'_',subj,'_',run,'.mat']);
save(savepath, 'ec_diff', '-v7.3');

gtm = 'eigenvector_cent';
param = [];
param.mask = gtm;
param.loc = 'max';
network_int_lcmv = vy_source_plot(ec_diff,template_mri,param,2);
savefig = fullfile(outputdir1,[mtd,'_1_',subj,'_',run]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

clear savepath
savepath{1} = fullfile(outputdir1,[mtd,'_2_',subj,'_',run]);
savepath{2} = fullfile(outputdir1,[mtd,'_3_',subj,'_',run]);
vy_mapvisualisation(network_int_lcmv,gtm,0.5, savepath,0);
% vy_mapvisualisation(network_int_lcmv,gtm,0.6, [],0);

savenii = fullfile(outputdir1,['n_',subj,'_',run,'.nii']);
vy_savenifti(network_int_lcmv, gtm, savenii);

%% parcellation - aal (132 rois)
% network_diff_lcmv.eigenvector_cent = (network_diff_lcmv.eigenvector_cent)./max(network_diff_lcmv.eigenvector_cent);
[~, data_intpar, coor] = vy_parcellate(ec_diff, atlas, gtm);
data_intpar.eigenvector_centdimord = 'chan';
% 
% savepath = fullfile(outputdir1,[mtd,'_par',subj,'_',run,'.mat']);
% save(savepath, 'data_intpar', '-v7.3');
% 
% clear savepath
% savepath{1} = fullfile(outputdir1,[mtd,'_par_1',subj,'_',run]);
% savepath{2} = fullfile(outputdir1,[mtd,'_par_2',subj,'_',run]);
% vy_mapvisualisation(data_intpar,gtm,0.3, savepath);
vy_mapvisualisation(data_intpar,gtm,0.6, []);
% 
% % save nii
% savenii = fullfile(outputdir1,['n_par',subj,'_',run,'.nii']);
% vy_savenifti(data_intpar, gtm, savenii);

%% Conn
conn_par.conn_thre = 0.7;
aedge =  cor_pst - cor_bsl;
ROI  = s_data2.pst.pos(s_data2.pst.inside,:);

aedge(isnan(aedge)) = 0;
tedge = (aedge.* double(aedge > conn_par.conn_thre.*max(aedge(:))));
figure, imagesc(tedge), colorbar

figure,
ft_plot_vol(individual_headmodel, 'facecolor', [0,0,0], 'edgecolor', 'none');
alpha 0.1;
view ([-196 56])
hold on
k = 1;
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
view([156,47]);
savepath = fullfile(outputdir1,['conn_temp',subj,'_',run]);
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
savepath = fullfile(outputdir1,[mtd,'_conn_',subj,'_',run,'.mat']);
save(savepath, 'data_intpar', 'individual_headmodel','-v7.3');

%% ROI summary
[ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, gtm);
disp(ROI_sel)
% savepath = fullfile(outputdir1,['n_ROIs_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
