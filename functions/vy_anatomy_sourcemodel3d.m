function sourcemodel = vy_anatomy_sourcemodel3d(mri, res)

if nargin<2
  res = 6;
end

cfg                 = [];
cfg.grid.warpmni    = 'yes';
cfg.grid.resolution = res;
cfg.grid.nonlinear  = 'yes';
cfg.mri             = mri;
sourcemodel         = ft_prepare_sourcemodel(cfg);
