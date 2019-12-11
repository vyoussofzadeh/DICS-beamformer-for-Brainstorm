function [freq,ff, psd,tapsmofrq] = vy_fft(cfg_mian, data)

%%
cfg              = [];
cfg.method       = 'mtmfft';
% cfg.output       = 'pow';
% cfg.pad          = 'nextpow2';
cfg.output       = 'fourier';
cfg.keeptrials   = 'yes';
cfg.foilim       = cfg_mian.foilim;
cfg.tapsmofrq    = cfg_mian.tapsmofrq;
cfg.taper        = cfg_mian.taper; %'hanning';
% cfg.taper        = 'dpss';
cfg.pad          = 4;
freq             = ft_freqanalysis(cfg, data);
psd = squeeze(mean(mean(abs(freq.fourierspctrm),2),1));
ff = linspace(1, cfg.foilim(2), length(psd));

if cfg_mian.plotflag ==1
    figure,plot(ff,psd)
    xlabel('Hz'); ylabel('psd')
%     max(psd)
end

tapsmofrq = cfg.tapsmofrq;

if cfg_mian.saveflag ==1
    hcp_write_figure([cfg_mian.savefile,'.png'], gcf, 'resolution', 300); 
end