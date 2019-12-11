function vy_ave_plot(cfg_main, data)


%% Global mean field power calculation for visualization purposes
cfg = [];
cfg.method = 'amplitude';
gmfp = ft_globalmeanfield(cfg, data);
%
figure;
pol = -1;     % correct polarity
scale = 10^6; % scale for eeg data micro volts
signal = scale*pol*data.avg; % add signle trials in a new value
% plot single trial together with global mean field power
h1 = plot(data.time,signal,'color',[0,0,0.5]);
hold on;
h2 = plot(data.time,scale*gmfp.avg,'color',[1,0,0],'linewidth',1);
%
legend([h1(1,1),h2],{'MEG','Gmean'});
grid on;
ylabel('MEG (fT)','Interpreter','Tex');
xlabel('Time (s)')
set(gca,'fontsize',18,'fontname','Century Gothic');

mx = max(max(signal));
mn = min(min(signal));
axis([data.time(1) data.time(end) mn mx])

% select time of interest for the source reconstruction later on
idx = find(data.time>0.6 & data.time<=0.9);
toi = data.time(idx);

[~,idxm] = max(max(abs(data.avg(:,idx))));
toi_mean_trial = toi(idxm);
%%
hcp_write_figure([cfg_main.savefile,'.png'], gcf, 'resolution', 300);

%%
cfg          = [];
cfg.fontsize = 6;
cfg.layout   = cfg_main.lay;
cfg.fontsize = 14;
% cfg.ylim     = [-5e-6 5e-6];
cfg.xlim     = [data.time(1) data.time(end)];
figure;
ft_multiplotER(cfg, data);

%%
% cfg            = [];
% cfg.zlim       = 'maxmin';
% cfg.comment    = 'xlim';
% cfg.commentpos = 'title';
% cfg.xlim       = [toi_mean_trial toi_mean_trial+0.01*toi_mean_trial];
% cfg.layout     = lay;
% cfg.fontsize   = 14;
%  
% figure; 
% ft_topoplotER(cfg, data);
