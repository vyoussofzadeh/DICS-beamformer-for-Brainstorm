function [n_data,fsb] = vy_notch(cfg_main, data)

if  ~isfield(cfg_main, 'fnotch') 
    cfg = [];
    cfg.savefile = [];
    cfg.saveflag = 2;
    cfg.foilim = [2 40];
    cfg.plotflag  = cfg_main.plotflag;
    cfg.tapsmofrq = 5;
    cfg.taper     = 'hanning';
    vy_fft(cfg, data);
    grid on
    grid minor
    title('Before band-stop filtering');
    fsb = input('Enter the sop-band frequency?');
else
    fsb = cfg_main.fnotch;
end
%
cfg = [];
cfg.bsfilter = 'yes';
%     cfg.bsfreq = [29 32]; % or whatever you deem appropriate
cfg.bsfreq = [fsb-1 fsb+1]; % or whatever you deem appropriate
%     cfg.bsfreq = [8 12;29 32]; % or whatever you deem appropriate
n_data = ft_preprocessing(cfg, data);
%     cfg.bsfreq = [2 12]; % or whatever you deem appropriate

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 40];
cfg.plotflag  = cfg_main.plotflag;
cfg.tapsmofrq = 5;
cfg.taper     = 'hanning';
vy_fft(cfg, n_data);
grid on
grid minor
title('After band-stop filtering');

end