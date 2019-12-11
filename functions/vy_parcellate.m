function [data_int, data_intpar, coor] = vy_parcellate(data, atlas, mask)

% atlas = ft_read_atlas(fullfile(atlas_path,'aal/ROI_MNI_V4.nii'));
% atlas = ft_read_atlas(fullfile(atlas_path,'brainnetome/BNA_MPM_thr25_1.25mm.nii'));


cfg = [];
cfg.parameter    = mask;
cfg.interpmethod = 'sph†ere_avg';
data_int  = ft_sourceinterpolate(cfg, data, atlas);
% data_int.dimord = 'pos';

cfg = [];
cfg.method      = 'mean';
% cfg.parcellation = 'parcellation';
data_intpar = ft_sourceparcellate(cfg, data_int, atlas);


% cfg = [];
% cfg.method = 'ortho';
% cfg.funparameter = 'pow';
% ft_sourceplot(cfg, data_intpar);
%
% cfg = [];
% cfg.funparameter = 'pow';
% cfg.funcolormap   = 'jet';
% cfg.method = 'surface';
% % cfg.projthresh     = 0.6;
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% ft_sourceplot(cfg, data_intpar)
%
% [L, idx] = sort(data_intpar.pow,'ascend');
% data_intpar.label(idx(1:10))'

%% MNI coordinates
coor = [];
pp = data_intpar.brainordinate;
for i=1:length(pp.tissuelabel)
    idx = pp.tissue == i;
    M = mean(pp.pos(idx,:),1);
    M = round(M*100)/100; % rounding into 2 decimals
    coor(i,:) = round(M);
end
% figure,
% plot3(coor(:,1),coor(:,2),coor(:,3),'b*');
% box off
% axis off
% axis tight
% set(gca,'color','none');

end

