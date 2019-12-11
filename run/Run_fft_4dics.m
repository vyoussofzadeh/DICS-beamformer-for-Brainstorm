cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 cfg_main.fmax];
cfg.plotflag  = 2;
cfg.tapsmofrq       = 1;
cfg.taper    = 'hanning';
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;

% % freq analysis - prepration for DICS source analysis
% f_data.bsl = vy_fft(ep_data.bsl, [2,40], 0,[],0); f_data.bsl.elec = sens;
% f_data.pst = vy_fft(ep_data.pst, [2,40], 0,[],0); f_data.pst.elec = sens;

% PSD - sensor space
psd_bsl = squeeze(mean(mean(abs(f_data.bsl.fourierspctrm),2),1));
psd_pst = squeeze(mean(mean(abs(f_data.pst.fourierspctrm),2),1));
ff = linspace(1, cfg_main.fmax, length(psd_pst));

% figure,plot(ff,psd_bsl)
% hold on
% % plot(ff,psd,'g')
% plot(ff,psd_pst,'r')
% hold on
% plot(ff,psd_pst - psd_bsl,'k')
% xlabel('Hz'); ylabel('psd'),legend({'bsl','pst','diff'})

% outputdir_dics = fullfile(cfg_main.outputdir,'dics');
outputdir_dics = cfg_main.outputdir;
if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end

% savepath = fullfile(outputdir_dics,['psd_',subj,'_',run,'.mat']);
% save(savepath, 'ff','psd_bsl','psd_pst', '-v7.3');
% savepath = fullfile(outputdir_dics,['fftcon_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

[a,b] = min(psd_pst - psd_bsl);
% f = ff(b); f = round(f);

f_sugg = round(cfg_main.freq_of_interest);
disp(['Suggested by TFR: ', num2str(f_sugg),'+-3Hz']);
disp('Select foi:')
%     if f_sugg > 12
%         f = f_sugg;
%     else
clear('input')
f = input('Freq of interest? '); 
tapsmofrq = 4;


% f = 22;     tapsmofrq = 4;
% f = 35;     tapsmofrq = 4;
% f=27;     tapsmofrq = 10;

%     end
%     f = f_sugg + tapsmofrq;
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [f f];
cfg.plotflag  = 2;
cfg.taper    = 'dpss'; cfg.tapsmofrq  = tapsmofrq;

if f < 4, cfg.tapsmofrq  = 1; cfg.taper    = 'hanning'; end
if toi(2,2)-toi(2,1) < 0.4, cfg.taper    = 'hanning'; end

[f_data.app,~,~,tapsmofrq] = vy_fft(cfg, ep_data.app); f_data.app.elec = cfg_main.sens;
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;

% f = 20;
% % Freq of interest - prepration for DICS source analysis
% f_data.app = vy_fft(ep_data.app, [f,f], 0,[],0); f_data.app.elec = sens;
% f_data.bsl = vy_fft(ep_data.bsl, [f,f], 0,[],0); f_data.bsl.elec = sens;
% f_data.pst = vy_fft(ep_data.pst, [f,f], 0,[],0); f_data.pst.elec = sens;