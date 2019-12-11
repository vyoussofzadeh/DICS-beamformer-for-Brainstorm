[s_data, s_data2] = vy_source(t_data, individual_grid, individual_headmodel);

%% reduce the source reconstructed data to the dominant orientation
cfg = [];
cfg.projectmom   = 'yes';
source_proj_post = ft_sourcedescriptives(cfg,s_data.pst);
source_proj_bsl  = ft_sourcedescriptives(cfg,s_data.bsl);

%% reduce memory demands and compute connectivity

% compute the sparse representation
source_sparse_post = ft_source2sparse(source_proj_post);
source_sparse_bsl  = ft_source2sparse(source_proj_bsl);

%%
data = source_sparse_bsl;
sizmom = size(data.avg.mom{data.inside(1)});
mom = zeros(size(data.pos,1), sizmom(2));
mom(data.inside, :) = cat(1, data.avg.mom{data.inside});
[nvox, nrpt]   = size(mom);
crsspctrm = (mom*mom')./nrpt;
% figure;imagesc(crsspctrm);
colorbar
mmom_bsl = mean(mom,1);

data = source_sparse_post;
sizmom = size(data.avg.mom{data.inside(1)});
mom = zeros(size(data.pos,1), sizmom(2));
mom(data.inside, :) = cat(1, data.avg.mom{data.inside});
[nvox, nrpt]   = size(mom);
crsspctrm = (mom*mom')./nrpt;
% figure;imagesc(crsspctrm);
colorbar
mmom_pst = mean(mom,1);

figure,
% plot(source_sparse_bsl.time,mmom_bsl)
plot(source_sparse_post.time,mmom_pst-mmom_bsl,'r')

restoredefaultpath
addpath(genpath(ft_old));
addpath(genpath([cd_org,'\functions']));
addpath(genpath([cd_org,'\data']));

% then compute connectivity
cfg = [];
cfg.method       = 'plv';
source_conn_post = ft_connectivityanalysis(cfg, source_proj_post);
source_conn_bsl  = ft_connectivityanalysis(cfg, source_proj_bsl);

source_conn_diff = source_conn_post;
source_conn_diff.plvspctrm = source_conn_post.plvspctrm - source_conn_bsl.plvspctrm ./ (source_conn_post.plvspctrm + source_conn_bsl.plvspctrm);
% source_conn_diff.plvspctrm = source_conn_post.plvspctrm;
% - source_conn_bsl.plvspctrm ./ (source_conn_post.plvspctrm + source_conn_bsl.plvspctrm);

%%
% figure;imagesc(source_conn_diff.plvspctrm);
%% network analysis
source_conn_full = (source_conn_diff);
source_conn_full.dimord    = 'pos_pos';

plv = source_conn_diff.plvspctrm(source_conn_diff.inside,source_conn_diff.inside);
% conn = tril(conn, -1);
figure;imagesc(plv);
% colormap(flipud(brewermap(64,'RdBu')));
grid off
box off
h = colorbar('eastoutside');
xlabel(h, 'PLV', 'FontSize', 14);


% figure;imagesc(source_conn_diff.cohspctrm);
% figure;imagesc(source_conn_diff.wpli_debiasedspctrm); 
% figure;imagesc(source_conn_diff.wplispctrm); 
colorbar
xlabel('voxels')
ylabel('voxels')
title('PLV')


% disp('=================');
% disp('degrees          = 1')
% disp('eigenvector_cent = 2')
% disp('betweenness      = 3');
% in2  = input('GT measure? ');
in2 = 2;
% computing graph metric
if in2 == 1
    gtm    = 'degrees';
elseif in2 == 2
    gtm    = 'eigenvector_cent';
elseif in2 == 3
    gtm    = 'betweenness';
elseif in2 == 4
    gtm    = 'efficiency_bin';
end

cfg = [];
cfg.method    = gtm;
% cfg.parameter = 'cohspctrm';
cfg.parameter = 'plvspctrm';
% cfg.parameter = 'wpli_debiasedspctrm';
% cfg.parameter = 'wplispctrm';
cfg.threshold = .5;
network = ft_networkanalysis(cfg,source_conn_full);

% if in2 ==2
%     network.eigenvector_cent(network.eigenvector_cent<0.0001) = 0;
% end

network.pos     = template_grid.pos;
network.dim     = template_grid.dim;
network.inside  = template_grid.inside;

if in2 ==1
        figure;bar(network.degrees(network.inside));
elseif in2 == 2
    figure;bar(network.eigenvector_cent(network.inside));
elseif in2 == 3
    figure;bar(network.betweenness(network.inside));
end
title('eigenvector centrality')

%%

% param = [];
% param.mask = gtm;
% param.loc = 'max';
% network_int = vy_source_plot(network,template_mri,param,2);
% 
% vy_mapvisualisation(network_int,gtm,0.6, [],0);

%% returning to the new ft!
addpath(genpath(fullfile(cd_org,'functions')));
vy_init(cd_org)

%%
network_diff_lcmv = network;
network_diff_lcmv.(gtm) = zscore(network_diff_lcmv.(gtm));

[mx,mx_idx] = max(network_diff_lcmv.eigenvector_cent(:));

outputdir1 = fullfile(outputdir,'network2');
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
vy_mapvisualisation(network_int_lcmv,gtm,0.6, savepath,0);
% vy_mapvisualisation(network_int_lcmv,gtm,0.4, [],0);

savenii = fullfile(outputdir1,['n_',subj,'_',run,'.nii']);
vy_savenifti(network_int_lcmv, gtm, savenii);

%%
ft_path = 'F:\My Matlab\My codes\My GitHub\fieldtrip_041718\fieldtrip-master';
% % % % ft_path = 'F:\My Matlab\My codes\GitHub\fieldtrip';
addpath(ft_path);
outputdir1 = fullfile(outputdir,'network2');
if exist(outputdir1, 'file') == 0
    mkdir(outputdir1);   %create the directory
end

network_int1 = network;
thre = 0.1;
network_int1.eigenvector_cent(network_int1.eigenvector_cent < thre*max(network_int1.eigenvector_cent(:)))=0;
% network_int1.eigenvector_cent(network_int1.eigenvector_cent > 0) = 1;
cfg = [];
cfg.filetype  = 'nifti';
cfg.parameter = gtm;
% cfg.filename  = '.\output';
cfg.filename  = fullfile(outputdir1,'output1.nii');
ft_volumewrite(cfg, network_int1)

%%
% vol  = individual_headmodel;
% aedge =  source_conn_diff.plvspctrm;
% tedge = zeros(size(aedge,1),size(aedge,2));
% tedge(mx_idx,:) = aedge(mx_idx,:);
% tedge(:,mx_idx) = aedge(:,mx_idx);

%%
vol  = individual_headmodel;
aedge =  source_conn_diff.plvspctrm;
tedge = (aedge.* double(aedge > 0.9.*max(max(aedge))));
for k = 1:length(tedge), lab{k} = num2str(k); end

Verticies = vol.bnd.pos;
ROI = source_conn_diff.pos(source_conn_diff.inside,:);
tedge = tedge(source_conn_diff.inside,source_conn_diff.inside);

figure,
plot3(Verticies(:,1),Verticies(:,2),Verticies(:,3),'color',[0.7,0.7,0.7]);
hold on
h = plot3(ROI(:,1),ROI(:,2),ROI(:,3),'.g');
% set(h(1),'MarkerEdgeColor',[1 0.48 0.30],'MarkerFaceColor','g')
set(h(1),'MarkerEdgeColor','g','MarkerFaceColor','g')
% set(h(2),'MarkerEdgeColor','none','MarkerFaceColor','g')
box off
set(gca,'color','none');
axis off
axis image
rotate3d on
hold on
view([-90,90])
% view ([180 90])

% figure,
% ft_plot_vol(vol, 'facecolor', 'none', 'edgecolor', [0.7,0.7,0.7]); alpha 0.5; camlight;
% hold on
% h = plot3(ROI(:,1),ROI(:,2),ROI(:,3),'.k');
% set(h(1),'MarkerEdgeColor',[1 0.48 0.30],'MarkerFaceColor','k')
% view ([-196 56])

if in2 == 1
    [L, Lidx] = sort(network.degrees(network.inside),'descend');
elseif in2 == 2
    [L, Lidx] = sort(network.eigenvector_cent(network.inside),'descend');
elseif in2 == 3
    [L, Lidx] = sort(network.betweenness(network.inside),'descend');
end

% tedge = verbs_wpli.cohspctrm;
nROI = 50;
for i = 1:nROI
    plot3(ROI(Lidx(i),1),ROI(Lidx(i),2),ROI(Lidx(i),3),'g.', 'MarkerSize',32);
end
for i = 1:length(tedge)
    for j = 1:length(tedge)
        if tedge(i,j)> 0.5
            p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
            p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
            pts = [p1; p2];
            line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
        end
    end
end

% cfg               = [];
% cfg.funcolormap   = 'jet';
% cfg.location      = 'max';
% cfg.funparameter  = gtm;
% cfg.method        = 'surface';
% cfg.surfinflated = vol.bnd;
% ft_sourceplot(cfg, network_int);
% alpha 0.4
% colorbar off
% % view([-90,90]);
% % view ([-196 56])
% view ([-263,20])

% vol2 = ft_convert_units(vol, 'mm');
figure,
ft_plot_vol(vol, 'facecolor', 'cortex', 'edgecolor', 'none'); alpha 0.5; camlight;
view ([-196 56])
% view ([27 39])

hold on
for i = 1:length(tedge)
    for j = 1:length(tedge)
        if tedge(i,j)> 0.5
            p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
            p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
            pts = [p1; p2];
            line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
        end
    end
end
hold on
for i = 1:nROI
    plot3(ROI(Lidx(i),1),ROI(Lidx(i),2),ROI(Lidx(i),3),'g.', 'MarkerSize',32);
    id(i) = Lidx(i);
end
set(gca,'color','none')
% title(['sum gtm for sub ', num2str(sub), ' is ', num2str(sum(L(1:nROI)))]);
% title(['']);

%%
restoredefaultpath
% close all
connpath = 'F:\My Matlab\Connectivity\Conn\conn';
addpath(connpath);
spm_path = 'F:\My Matlab\SPM\spm12_4\spm12\spm12';
addpath(genpath(spm_path))
h = get(0, 'Children');
if isempty(findobj(h,'tag','CONN functional connectivity toolbox'))
    conn
end

ROI = network.pos(network.inside,:);

% ROI2 = cor2mni(ROI, template_mri.transform );
%  ROI2 = cor2mni(ROI, network_int.initial );
conn_mesh_display('', '', '', ROI(id,:), tedge(id,id), .1);

%%
% outputdir1 = fullfile(outputdir,'network2');
% filenameSURF = fullfile('output1.nii');
% filenameVOL = [];
% FSfolder = fullfile(connpath,'\utils\surf');
% sphplots = [];
% connplots = [];
% facealpha = 1;
% position = [-1 0  0];
% conn_mesh_display(filenameSURF,[],FSfolder);

