function vy_savenifti(input,funparameter,savename)

cfg = [];
cfg.filetype  = 'nifti';
cfg.datatype   = 'uint8'; %'float';
cfg.parameter = funparameter;
cfg.filename  = savename;
ft_volumewrite(cfg, input)

disp(['nii was saved as, ',savename])
