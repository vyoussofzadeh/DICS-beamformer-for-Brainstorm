function vy_mapvisualisation(cfg_main, input)

% ,funparameter,thre,savepath

cfg                = [];
cfg.method         = 'ortho';
cfg.funparameter   = cfg_main.mask;
cfg.funcolorlim    = 'maxabs';
cfg.opacitymap     = 'rampup';
cfg.crosshair      = 'no';
cfg.camlight       = 'no';
cfg.funcolormap    =  brewermap(256, '*RdYlBu');
cfg.projthresh     = cfg_main.thre;

cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.surfinflated   = 'surface_inflated_right.mat';
ft_sourceplot(cfg, input);
% view ([-70 20 50])
% light ('Position',[-70 20 50]);
view([90 0]);
camlight; material dull;
set(gcf,'name',cfg_main.subj,'numbertitle','off')
% title('test');

if ~isempty(cfg_main.savepath)==1
    %     hcp_write_figure([savepath{1},'.png'], gcf, 'resolution', 300);
    pause(0.5)
%     print(cfg_main.savepath{1},'-depsc')
    print(cfg_main.savepath{1},'-dpng')    
end

ft_sourceplot(cfg, input);
% view ([70 20 50])
% light ('Position',[70 20 50])
view([-90 0]); camlight; material dull;
set(gcf,'name',cfg_main.subj,'numbertitle','off')

if ~isempty(cfg_main.savepath)==1
    %     hcp_write_figure([savepath{2},'.png'], gcf, 'resolution', 300);
    pause(0.5)
%     print(cfg_main.savepath{2},'-depsc')
    print(cfg_main.savepath{2},'-dpng')
end



