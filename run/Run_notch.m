cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 40];
cfg.plotflag  = 1;
cfg.tapsmofrq = 5;
cfg.taper     = 'hanning';
[freq,ff,psd] = vy_fft(cfg, cln_data);
grid on
grid minor
title(['Before band-stop filtering-',subj]);

%- finding notch freq
idx = find(freq.freq ==30); TF = islocalmax(psd); TF(1:idx-5) = 0;
hold on
plot(ff(TF),psd(TF),'r*')

idx2 = find(TF == 1);
if length(idx2)==1
    fsb = round(ff(idx2));
else
    [val, idx] = max(psd(TF));
    fff = ff(TF);
    %             fsb = input(['Enter the sop-band frequency for ', subj,'?']);
    fsb = round(fff(idx));
end
disp([num2str(fsb),'Hz freq was selected for notch filteting']);

%%
cfg = [];
cfg.bsfilter = 'yes';
%     cfg.bsfreq = [29 32]; % or whatever you deem appropriate
cfg.bsfreq = [fsb-1 fsb+1]; % or whatever you deem appropriate
%     cfg.bsfreq = [8 12;29 32]; % or whatever you deem appropriate
cln_data = ft_preprocessing(cfg, cln_data);
%     cfg.bsfreq = [2 12]; % or whatever you deem appropriate

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 40];
cfg.plotflag  = 1;
cfg.tapsmofrq = 5;
cfg.taper     = 'hanning';
vy_fft(cfg, cln_data);
grid on
grid minor
title(['After band-stop filtering-',subj]);