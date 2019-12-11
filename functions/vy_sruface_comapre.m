function vy_sruface_comapre(meg, meg_thresh, fmri, fmri_thresh)

%- meg
s_vol = ft_read_mri(niifile1);
s_vol_thre = vy_vol_thresh(s_vol,projthresh);

cfg = [];
cfg.filetype  = 'nifti';
cfg.datatype   = 'uint8'; %'float';
cfg.parameter = 'anatomy';
cfg.filename  = vol_name;
ft_volumewrite(cfg, s_vol_thre)

conn_mesh_display(vol_name, '');



end