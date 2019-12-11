function sourceint_pow = vy_parcellate_plot(data_intpar,coor, name)

%- surface visualisation
% cfg               = [];
% cfg.method        = 'surface';
% cfg.funparameter  = 'anatomy';
% cfg.maskparameter = cfg.funparameter;
% cfg.funcolorlim   = [-0.3 0.3];
% cfg.opacitylim    = [-8 8];
% cfg.opacitymap    = 'rampup';
% cfg.funcolormap   = 'jet';
% cfg.colorbar      = 'yes';
% cfg.projthresh     = 0.6;
% ft_sourceplot(cfg, data_intpar);
% view ([-100 0 0]), light ('Position',[-100 0 0])

sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
cfg1 = [];
cfg1.sourceunits  = 'mm';
cfg1.parameter    = 'anatomy';
cfg1.downsample   = 1;
% cfg.interpmethod = 'sphere_avg';
% cfg1.coordsys     = 'mni';
sourceint_pow = ft_sourceinterpolate(cfg1, data_intpar, ft_read_mri(sMRI, 'format', 'nifti_spm'));

sourceint_pow.anatomy(isnan(sourceint_pow.anatomy(:))) = 0;

clear savepath
savepath{1} = [name,'_left'];
savepath{2} = [name,'_right'];
cfg = [];
cfg.subj = 'par';
cfg.mask = 'anatomy';
cfg.thre = 0.7;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, sourceint_pow);
% vy_mapvisualisation(sourceint_pow,'anatomy',0.7, savepath);

vy_ROI_report(data_intpar,.7, coor, 'anatomy');
savepath = ['ROI_',name];
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300)

