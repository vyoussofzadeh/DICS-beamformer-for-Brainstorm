speechonset_time = tt(idx_good);
% speechonset_time = speechonset_time(ia);
L = length(speechonset_time);

speechonset_sample = ipoints_all(idx_good);
% speechonset_sample = speechonset_sample(ia);

figure,
plot(speechonset_time, '*'),
hline(thre,'b',['mean:', num2str(thre),'sec']),
box off;
set(gca,'color','none');
title('Speech onset'),
ylabel('Time (sec)');
xlabel('Trials');
set(gca,'Xtick', 1:L,'XtickLabel',1:L);
xlim([1 L]);
set(gca,'FontSize',10,'XTickLabelRotation',90);


%%
ep_data.all = cln_data;
ep_data.pst = [];

%% PST
tinterval = 700;
trl_res = [];
for i = 1:length(cln_data.trial)   
    trl = cln_data.trial{i};
    trl_res.trl{i} = trl(:, speechonset_sample(i)-tinterval:speechonset_sample(i));
    ttrl = cln_data.time{i};
%     trl_res.time{i} = ttrl(:, speechonset_sample(1)-1000:speechonset_sample(1));
    trl_res.time{i} = linspace(-tinterval/1000,0,tinterval+1);
end

ep_data.pst = cln_data;
ep_data.pst.trial = trl_res.trl;
ep_data.pst.time = trl_res.time;

%% BSL
tinterval = 300;
trl_res = [];
for i = 1:length(cln_data.trial)
    trl = cln_data.trial{i};
    idx = find(cln_data.time{1} == 0);     idx1 = find(cln_data.time{1} == -tinterval/1000);
    trl_res.trl{i} = trl(:, idx1:idx);
    ttrl = cln_data.time{i};
    %     trl_res.time{i} = ttrl(:, speechonset_sample(1)-1000:speechonset_sample(1));
    trl_res.time{i} = linspace(-tinterval/1000,0,tinterval+1);
end

ep_data.bsl = cln_data;
ep_data.bsl.trial = trl_res.trl;
ep_data.bsl.time = trl_res.time;


%- Appending data
cfg = [];
ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);

%% Data-covarinace estimation
t_data = vy_timelock(ep_data);

%%
toi = [-700, 0; -300, 0];

%% tfr analysis
tinterval = 1000;
trl_res = [];
for i = 1:length(cln_data.trial)   
    trl = cln_data.trial{i};
    trl_res.trl{i} = trl(:, speechonset_sample(i)-tinterval:speechonset_sample(i));
    ttrl = cln_data.time{i};
%     trl_res.time{i} = ttrl(:, speechonset_sample(1)-1000:speechonset_sample(1));
    trl_res.time{i} = linspace(-tinterval/1000,0,tinterval+1);
end

ep_test = cln_data;
ep_test.trial = trl_res.trl;
ep_test.time = trl_res.time;


toi = [ ep_test.time{1}(1),  ep_test.time{1}(end)];

cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
% cfg.taper      = 'dpss';
cfg.foi        = 1:3:40;
cfg.keeptrials = 'yes';
cfg.t_ftimwin  = 3./cfg.foi;
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.tapsmofrq  = 0.8 *cfg.foi;
cfg.toi        =  toi(1):0.05:toi(2);
tfr        = ft_freqanalysis(cfg, ep_test);


% toi_bsl = [-0.3,0];
% cfg = [];
% cfg.output     = 'pow';
% cfg.channel    = 'all';
% cfg.method     = 'mtmconvol';
% cfg.taper      = 'hanning';
% % cfg.taper      = 'dpss';
% cfg.foi        = 1:3:40;
% cfg.keeptrials = 'yes';
% cfg.t_ftimwin  = 3./cfg.foi;
% % cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
% cfg.tapsmofrq  = 0.8 *cfg.foi;
% cfg.toi        =  toi_bsl(1):0.05:toi_bsl(2);
% tfr_bsl        = ft_freqanalysis(cfg, ep_data.bsl);

cfg = [];
freq_avg = ft_freqdescriptives(cfg, tfr);
% freq_avg_bsl = ft_freqdescriptives(cfg, tfr_bsl);

cfg = [];
cfg.baseline = [-0.3 0];
cfg.baselinetype = 'db'; % Use decibel contrast here
freq_avg_bsl = ft_freqbaseline(cfg, freq_avg);

meanpow = squeeze(mean(freq_avg_bsl.powspctrm, 1));
% freq_avg_bsl.powspctrm(freq_avg_bsl.powspctrm > 0) = 0; 

tim_interp = linspace(toi(1), toi(2), 512);
freq_interp = linspace(1, 40, 512);

% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(tfr.time, tfr.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

meanpow(isnan(meanpow))=0;


% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, tim_grid_interp, freq_grid_interp, 'spline');

%
pow_interp1  = pow_interp(50:end,50:end);
tim_interp1  = tim_interp(50:end);
freq_interp1 = freq_interp(50:end);
% pow_interp(pow_interp>0)=0;

%
[~,idx] = min(pow_interp1(:));
[row,col] = ind2sub(size(pow_interp1),idx);

time_of_interest = tim_interp1(col);
freq_of_interest = freq_interp1(row);

% [a, idx] = min(pow_at_toi(20:end,:));
% freq_interp(idx)

% time_of_interest = 1.6;
% freq_of_interest = 20;
timind = nearest(tim_interp, time_of_interest);
freqind = nearest(freq_interp, freq_of_interest);
pow_at_toi = pow_interp(:,timind);
pow_at_foi = pow_interp(freqind,:);

%
figure();
ax_main  = axes('Position', [0.1 0.2 0.55 0.55]);
ax_right = axes('Position', [0.7 0.2 0.1 0.55]);
ax_top   = axes('Position', [0.1 0.8 0.55 0.1]);

%
axes(ax_main);
im_main = imagesc(tim_interp, freq_interp, pow_interp);
% note we're storing a handle to the image im_main, needed later on
xlim([toi(1), toi(2)]);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
clim = max(abs(meanpow(:)));
caxis([-clim clim]);
colormap(brewermap(256, '*RdYlBu'));
hold on;
plot(zeros(size(freq_interp)), freq_interp, 'k:');

%
axes(ax_top);
area(tim_interp, pow_at_foi,...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
xlim([toi(1), toi(2)]);
ylim([-clim clim]);
box off;
ax_top.XTickLabel = [];
ylabel('Power (dB)');
hold on;
plot([0 0], [-clim clim], 'k:');


axes(ax_right);
area(freq_interp, pow_at_toi,...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
view([270 90]); % this rotates the plot
ax_right.YDir = 'reverse';
ylim([-clim clim]);
box off;
ax_right.XTickLabel = [];
ylabel('Power (dB)');

%
h = colorbar(ax_main, 'manual', 'Position', [0.85 0.2 0.05 0.55]);
ylabel(h, 'Power vs baseline (dB)');


%
% Main plot:
axes(ax_main);
plot(ones(size(freq_interp))*time_of_interest, freq_interp,...
    'Color', [0 0 0 0.1], 'LineWidth', 3);
plot(tim_interp, ones(size(tim_interp))*freq_of_interest,...
    'Color', [0 0 0 0.1], 'LineWidth', 3);

% Marginals:
axes(ax_top);
plot([time_of_interest time_of_interest], [0 clim],...
    'Color', [0 0 0 0.1], 'LineWidth', 3);
axes(ax_right);
hold on;
plot([freq_of_interest freq_of_interest], [0 clim],...
    'Color', [0 0 0 0.1], 'LineWidth', 3);

% set(gca, 'Xcolor', 'k', 'Ycolor', 'k')
% set(gca, 'YAxisLocation', 'right')


%%
% cfg = [];
% cfg.baseline = [-0.3 0];
% cfg.baselinetype = 'absolute';
% freq_bsl = ft_freqbaseline(cfg, tfr);
% 
% cfg = [];
% cfg.variance = 'yes';
% freq_sem = ft_freqdescriptives(cfg, freq_bsl);
% 
% tscore = freq_sem.powspctrm./ freq_sem.powspctrmsem;
% 
% % Average the t-score over our channels:
% tscore = squeeze(mean(tscore, 1));
% 
% % figure();
% % imagesc(tfr.time, tfr.freq, tscore);
% % axis xy;
% % colorbar();
% 
% tscore(isnan(tscore))=0;
% 
% tscore_interp = interp2(tim_grid_orig, freq_grid_orig, tscore,...
%     tim_grid_interp, freq_grid_interp, 'spline');
% 
% alpha = 0.01;
% tcrit = tinv(1-alpha/2, size(tfr.powspctrm, 1)-1);
% 
% opacity = abs(tscore_interp) / tcrit;
% opacity(opacity > 1) = 1;
% 
% im_main.AlphaData = opacity;
% 
% if ~isempty(cfg_main.savefile)
% hcp_write_figure([cfg_main.savefile,'.png'], gcf, 'resolution', 300);
% end

%%
% cfg = [];
% cfg.savepath = fullfile(outd.sub,'Freq');
% cfg.savefile = fullfile(savepath,[stag,'_tfr2_',subj]);
% cfg.toi = [tfr.time(1), tfr.time(end)];
% cfg.fmax = fmax;
% [time_of_interest,freq_of_interest] = vy_tfr_plot(cfg, tfr);
