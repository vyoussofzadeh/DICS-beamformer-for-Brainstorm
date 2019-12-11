function [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level)


[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(level, info.chansel, info.trlsel, strcmp(info.metric, 'min'));

if info.pflag == 1
    figure,
    axis ij;
    ymax = max(maxperchan); ymin = min(maxperchan); xmax = info.nchan;
    subplot 121,
    plot(maxperchan,'.'), xlabel('Channel number'), ylabel(info.metric),
    axis([0.5 xmax+0.5 0.8*ymin 1.2*ymax]);
    title('Channel')
    subplot 122,
    plot(maxpertrl,'.'), xlabel('Trial number'), ylabel(info.metric);
    xmax = info.ntrl; ymax = max(maxpertrl); ymin = min(maxpertrl);
    axis([0.5 xmax+0.5 (1-sign(ymin)*0.2)*ymin (1+sign(ymax)*0.2)*ymax]);
    title('Trial')
end


function [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(level, chansel, trlsel, minflag)
if minflag
    % take the negative maximum, i.e. the minimum
    level = -1 * level;
end
% determine the maximum value
maxperchan_all = max(level, [], 2);
maxpertrl_all  = max(level, [], 1);
% determine the maximum value over the remaining selection
level(~chansel, :) = nan;
level(:, ~trlsel)  = nan;
maxperchan     = max(level, [], 2);
maxpertrl      = max(level, [], 1);
if minflag
    maxperchan     = -1 * maxperchan;
    maxpertrl      = -1 * maxpertrl;
    maxperchan_all = -1 * maxperchan_all;
    maxpertrl_all  = -1 * maxpertrl_all;
    level          = -1 * level;
end
