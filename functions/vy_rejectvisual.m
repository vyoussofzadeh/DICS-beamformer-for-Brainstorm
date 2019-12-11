function data_fix = vy_rejectvisual(data_in, lay, savepath,saveflag)

if exist(savepath, 'file') == 2
    load(savepath)
else
    cfg = [];
    cfg.metric = 'zvalue';  % use by default zvalue method
    cfg.latency = [-400,900];
    cfg.layout   = lay;   % this allows for plotting individual trials
    data_fix   = ft_rejectvisual(cfg, data_in);
    if saveflag ==1
        save(savepath, 'data_fix', '-v7.3');
    end
end