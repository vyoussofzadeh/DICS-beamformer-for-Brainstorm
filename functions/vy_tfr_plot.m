function [time_of_interest,freq_of_interest] = vy_tfr_plot(cfg_main, tfr)

%%
% First compute the average over trials:
cfg = [];
freq_avg = ft_freqdescriptives(cfg, tfr);

% And baseline-correct the average:
cfg = [];
cfg.baseline = [-0.3 0];
cfg.baselinetype = 'db'; % Use decibel contrast here
freq_avg_bsl = ft_freqbaseline(cfg, freq_avg);

freq_avg_bsl.powspctrm(isnan(freq_avg_bsl.powspctrm))=0;
meanpow = squeeze(mean(freq_avg_bsl.powspctrm, 1));

%     figure();
%     imagesc(tfr.time, tfr.freq, meanpow);
%     xlim([-0.5 2]); % to hide the NaN values at the edges of the TFR
%     axis xy; % to have the y-axis extend upwards
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
%     clim = max(abs(meanpow(:)));
%     caxis([-clim clim]);
% 
% % addpath('thirdparty\brewermap'); % or wherever you put brewermap
% 
%     colormap(brewermap(256, '*RdYlBu'));
%     h = colorbar();
%     ylabel(h, 'Power vs baseline (dB)');

%%
% The finer time and frequency axes:
% tim_interp = linspace(-0.5, 2, 512);
% freq_interp = linspace(1, 40, 512);

tim_interp = linspace(cfg_main.toi(1), cfg_main.toi(2), 512);
freq_interp = linspace(1, cfg_main.fmax, 512);

% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(tfr.time, tfr.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% meanpow(isnan(meanpow))=0;


% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, tim_grid_interp, freq_grid_interp, 'spline');

%%
% figure();
% imagesc(tim_interp, freq_interp, pow_interp);
% xlim([-0.5 2]);
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% clim = max(abs(meanpow(:)));
% caxis([-clim clim]);
% colormap(brewermap(256, '*RdYlBu'));
% h = colorbar();
% ylabel(h, 'Power vs baseline (dB)');
% 
% % Let's also add a line at t = 0s for clarity:
% hold on;
% plot(zeros(size(freq_interp)), freq_interp, 'k:');

%%
% while n==1
pow_interp1  = pow_interp(50:end,50:end);
tim_interp1  = tim_interp(50:end);
freq_interp1 = freq_interp(50:end);

%%
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

%%
figure();
ax_main  = axes('Position', [0.1 0.2 0.55 0.55]);
ax_right = axes('Position', [0.7 0.2 0.1 0.55]);
ax_top   = axes('Position', [0.1 0.8 0.55 0.1]);

%%
axes(ax_main);
im_main = imagesc(tim_interp, freq_interp, pow_interp);
% note we're storing a handle to the image im_main, needed later on
xlim([cfg_main.toi(1), cfg_main.toi(2)]);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
clim = max(abs(meanpow(:)));
caxis([-clim clim]);
colormap(brewermap(256, '*RdYlBu'));
hold on;
plot(zeros(size(freq_interp)), freq_interp, 'k:');

%%
%     min_power = min(pow_at_foi);
%     idx_time = find(pow_at_foi == min_power);
%     time_of_interest = tim_interp(idx_time);
%
%     min_power = min(pow_at_toi);
%     idx_freq = find(pow_at_toi == min_power);
%     freq_of_interest = freq_interp(idx_freq);

%%
axes(ax_top);
area(tim_interp, pow_at_foi,...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
xlim([cfg_main.toi(1), cfg_main.toi(2)]);
ylim([-clim clim]);
box off;
ax_top.XTickLabel = [];
ylabel('Power (dB)');
hold on;
plot([0 0], [-clim clim], 'k:');

%%
% [a, idx] = min(pow_at_toi(20:end,:));
% freq_interp(idx)
% % 
% figure,plot(freq_interp,pow_at_toi)

%%
axes(ax_right);
area(freq_interp, pow_at_toi,...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
view([270 90]); % this rotates the plot
ax_right.YDir = 'reverse';
ylim([-clim clim]);
box off;
ax_right.XTickLabel = [];
ylabel('Power (dB)');

%%
h = colorbar(ax_main, 'manual', 'Position', [0.85 0.2 0.05 0.55]);
ylabel(h, 'Power vs baseline (dB)');

%%
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
cfg = [];
cfg.baseline = [-0.3 0];
cfg.baselinetype = 'absolute';
freq_bsl = ft_freqbaseline(cfg, tfr);

cfg = [];
cfg.variance = 'yes';
freq_sem = ft_freqdescriptives(cfg, freq_bsl);

tscore = freq_sem.powspctrm./ freq_sem.powspctrmsem;

% Average the t-score over our channels:
tscore = squeeze(mean(tscore, 1));

% figure();
% imagesc(tfr.time, tfr.freq, tscore);
% axis xy;
% colorbar();

tscore(isnan(tscore))=0;

tscore_interp = interp2(tim_grid_orig, freq_grid_orig, tscore,...
    tim_grid_interp, freq_grid_interp, 'spline');

alpha = 0.01;
tcrit = tinv(1-alpha/2, size(tfr.powspctrm, 1)-1);

opacity = abs(tscore_interp) / tcrit;
opacity(opacity > 1) = 1;

% im_main.AlphaData = opacity;

if ~isempty(cfg_main.savefile)
% hcp_write_figure([cfg_main.savefile,'.png'], gcf, 'resolution', 300);
saveas(gcf,[cfg_main.savefile,'.png'])
end

%%
% [~,idx] = min(tscore_interp(:));
% [row,col] = ind2sub(size(tscore_interp),idx);
% 
% time_of_interest = tim_interp(col);
% freq_of_interest = freq_interp(row);
% % time_of_interest = 1.6;
% % freq_of_interest = 20;
% timind = nearest(tim_interp, time_of_interest);
% freqind = nearest(freq_interp, freq_of_interest);
% t_at_toi = tscore_interp(:,timind);
% t_at_foi = tscore_interp(freqind,:);

%%
% figure();
% ax_main = axes('Position', [0.1 0.2 0.55 0.55]);
% ax_right = axes('Position', [0.7 0.2 0.1 0.55]);
% ax_top = axes('Position', [0.1 0.8 0.55 0.1]);
% 
% axes(ax_main);
% im_main = imagesc(tim_interp, freq_interp, tscore_interp);
% % note we're storing a handle to the image im_main, needed later on
% xlim([-0.5 2]);
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% clim = max(abs(tscore(:)));
% caxis([-clim clim]);
% colormap(brewermap(256, '*RdYlBu'));
% hold on;
% plot(zeros(size(freq_interp)), freq_interp, 'k:');
% 
% axes(ax_top);
% area(tim_interp, t_at_foi,...
%     'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
% xlim([-0.5 2]);
% ylim([-clim clim]);
% box off;
% ax_top.XTickLabel = [];
% ylabel('Power (dB)');
% hold on;
% plot([0 0], [-clim clim], 'k:');
% 
% axes(ax_right);
% area(freq_interp, t_at_toi,...
%     'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
% view([270 90]); % this rotates the plot
% ax_right.YDir = 'reverse';
% ylim([-clim clim]);
% box off;
% ax_right.XTickLabel = [];
% ylabel('Power (dB)');
% 
% h = colorbar(ax_main, 'manual', 'Position', [0.85 0.2 0.05 0.55]);
% ylabel(h, 't-value vs baseline (dB)');

end
