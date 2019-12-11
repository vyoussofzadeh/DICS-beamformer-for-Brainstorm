function vy_surfce_vis(s_vol,vol_name, opt)

if opt.savenii == 1
    cfg = [];
    cfg.filetype  = 'nifti';
    cfg.datatype   = 'uint8'; %'float';
    cfg.parameter = 'anatomy';
    cfg.filename  = vol_name;
    ft_volumewrite(cfg, s_vol);
end

conn_mesh_display(vol_name, '');
title(['\fontsize{16} {\color{magenta}',[opt.tsk,', ',opt.subj,'-', opt.run],'}'],'interpreter','tex');
pause(1.5);

savefig = fullfile(opt.savedir,[opt.tsk, '_', opt.subj,'_',opt.run,'.jpg']);

conn_print(savefig,...
    '-nogui',...
    '-mosaic',...
    get(findobj(gcf,'label','Left view'),'callback'),...
    get(findobj(gcf,'label','Left medial view'),'callback'),...
    get(findobj(gcf,'label','Right view'),'callback'),...
    get(findobj(gcf,'label','Right medial view'),'callback'));

end