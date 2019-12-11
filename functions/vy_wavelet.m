function w_data = vy_wavelet(data, lay, savepath,foi)


cfg = [];
cfg.method = 'wavelet'; %mtmconvol
cfg.output = 'powandcsd';
cfg.foi = foi;
cfg.toi = -0.400:0.020:1.000;
w_data = ft_freqanalysis(cfg, data);
save(savepath, 'w_data', '-v7.3');

%% Plotting wavelet
figure,
cfg = [];
cfg.layout = lay;
cfg.baseline = [-inf 0];
cfg.baselinetype = 'relative';
ft_multiplotTFR(cfg, w_data);

if isempty(savepath) == 0
    hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300); 
end
