function vy_surfmap(cfg_main, input)

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = cfg_main.maskparam;
cfg.funcolorlim = 'maxabs';
cfg.opacitymap = 'rampup';
cfg.crosshair = 'no';
cfg.camlight       = 'no';
cfg.funcolormap =  brewermap(256, '*RdYlBu');
cfg.projthresh     = cfg_main.projthresh;

cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.surfinflated   = 'surface_inflated_right.mat';
ft_sourceplot(cfg, input);
% view ([-70 20 50])
% light ('Position',[-70 20 50]);
view([90 0]);
camlight; material dull;
% title('test');

print -dpng

if ~isempty(cfg_main.save.savepath)==1
%     hcp_write_figure([savepath{1},'.png'], gcf, 'resolution', 300);
    export_fig(cfg_main.save.savepath{1}, cfg_main.save.saveformat, cfg_main.save.pixdim)
end

ft_sourceplot(cfg, input);
% view ([70 20 50])
% light ('Position',[70 20 50])
view([-90 0]); camlight; material dull;

if ~isempty(cfg_main.save.savepath)==1
%     hcp_write_figure([savepath{2},'.png'], gcf, 'resolution', 300);
    export_fig(cfg_main.save.savepath{2}, cfg_main.save.saveformat, cfg_main.save.pixdim)
end