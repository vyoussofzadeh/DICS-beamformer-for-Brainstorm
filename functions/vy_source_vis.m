function vy_source_vis(funparameter,input,tit,savename,thre, saveROI,savepath)


%%
% r = [1 0 0];       %# start
% w = [1 1 1];    %# middle
% b = [0 0 1];       %# end
%
% %# colormap of size 64-by-3, ranging from red -> white -> blue
% k = 90;
% c1 = zeros(k,3); c2 = zeros(k,3);
% for i=1:3
%     c1(:,i) = linspace(b(i), w(i), k);
%     c2(:,i) = linspace(w(i), r(i), k);
% end
% c = [c1;c2];
%%

% r = [1 0 0];       %# start
% % w = [1 1 1];    %# middle
% blk = [0 0 0];       %# end
% w = [1 1 1];       %# end
%
% k = 256;
% c = zeros(k,3);
% for i=1:3
% %     c(:,i) = linspace(blk(i), r(i), k);
%     c(:,i) = linspace(w(i), r(i), k);
% %     c2(:,i) = linspace(w(i), r(i), k);
% end
% % c = [c1;c2];
%
% colormap(c), colorbar
%%
% c = jet(256);

cfg               = [];
cfg.method        = 'ortho';
% cfg.funcolormap   = c;
% cfg.location      = 'max';
cfg.funparameter  = funparameter;
cfg.funcolorlim = 'maxabs';
cfg.opacitymap = 'rampup';
cfg.crosshair = 'no';
cfg.camlight       = 'no';
cfg.projthresh     = thre;
% ft_sourceplot(cfg, input);
% set(gca,'color','none');
% title(tit)


% ft_plot_ortho(input)

cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.funcolormap = flipud(brewermap(64,'RdBu'));
% cfg.funcolormap   = jet;
figure,
subplot 221
ft_sourceplot(cfg, input);
view ([-70 20 50])
% view ([-90 0])
% view ([70 20 50])
% light ('Position',[-70 20 50])
% colormap jet
% title(tit)
hcp_write_figure([savepath{1},'.png'], gcf, 'resolution', 300);


% ft_sourceplot(cfg, input);
% view ([-70 20 50])
% % view ([-90 0])
% % view ([70 20 50])
% light ('Position',[-70 20 50])
%%
% if isempty(savename)== 0
%     path = '.\output';
%     saveas(gca, fullfile(path, [savename,'_L']), 'jpeg');
% end
% input.eigenvector_cent (abs(input.eigenvector_cent) < thre*max(abs(input.eigenvector_cent))) = 0;

ft_sourceplot(cfg, input);
view ([70 20 50])
% view ([90 0])
light ('Position',[70 20 50])
hcp_write_figure([savepath{2},'.png'], gcf, 'resolution', 300);

if isempty(savename) == 0
    %     saveas(gca, fullfile(path, [savename,'_R']), 'jpeg');
    cfg = [];
    cfg.filetype  = 'nifti';
    cfg.parameter = funparameter;
    cfg.filename  = savename;
    ft_volumewrite(cfg, input)
end
%% Visualise and export first 5 ROIs as nii
% addpath(genpath('E:\My Matlab\My codes\My GitHub\My tool\MEG_ft_pipline\Time domain\External'));
% addpath(genpath('E:\My Matlab\My codes\My GitHub\fieldtrip_new\fieldtrip-20170321'));
% warning off
if saveROI ==1
    input1 = input;
    idx = find(input1.eigenvector_cent> 0);
    [idx2,b] = sort(input1.eigenvector_cent(idx),'descend');
    for i=1:length(idx)
        ll = input1.label(idx(b(i)))
        r = zeros(116,1);
        r(idx(b(i))) = 1;
        input1.eigenvector_cent = r;
        %     ft_mapvisualisation('eigenvector_cent',mm1, [],'all_parc');
        
        cfg = [];
        cfg.filetype  = 'nifti';
        cfg.parameter = 'eigenvector_cent';
        cfg.filename  = [savename,'_',num2str(i),'_',ll{1}];
        ft_volumewrite(cfg, input1)
        
    end
end


