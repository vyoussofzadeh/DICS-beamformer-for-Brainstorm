function [ri,outliers] = vy_outlier_baseline(dat, baseline)

%---------------------------------------------------------------
% This will work.
% Compute the median absolute difference
% baseline_dat = dat(baseline,:);

dat(isnan(dat))=0;
meanValue = mean(dat,1);
% 
for i=1:size(dat,1)
    r(i) = corr2(dat(i,:),meanValue);
end
[ri, outliers] = sort(r);

% kurt = kurtosis(dat, [], 2);
% [ri, outliers] = sort(kurt);
% outliers = outliers';
% ri = ri';


label = {'var', 'corr with mean', 'max', 'maxabs', 'range','kurtosis','1/var'};
% 'var'
level(:, 1) = nanstd(dat, [], 2).^2; 
% 'corr'
level(:, 2) = r;
% 'max'
level(:, 3) = nanmax(dat, [], 2);
% 'maxabs'
level(:, 4) = nanmax(abs(dat), [], 2);
% 'range'
level(:, 5) = nanmax(dat, [], 2) - nanmin(dat, [], 2);
% 'kurtosis'
level(:, 6) = kurtosis(dat, [], 2);
% '1/var'
level(:, 7) = 1./(nanstd(dat, [], 2).^2);

figure,
for i=1:7
    subplot(4,2,i)
    plot(level(:,i),'.');title(label{i})
end

% figure,
% plot(r,'o', 'MarkerFaceColor', 'b');
% title('corr')
% grid on
% box off
% set(gcf, 'Color', 'None')
% warning('off','MATLAB:hg:ColorSpec_None')


% % Compute the absolute differences.  It will be a vector.
% absoluteDeviation = abs(data - meanValue);
% % Compute the median of the absolute differences
% mad = median(absoluteDeviation,1);
% % Find outliers.  They're outliers if the absolute difference
% % is more than some factor times the mad value.
% sensitivityFactor = 6; % Whatever you want.
% thresholdValue = sensitivityFactor * mad;
% outlierIndexes = abs(absoluteDeviation) > thresholdValue;
% % Extract outlier values:
% outliers = double(outlierIndexes);
% % Extract non-outlier values:
% nonOutliers = double(~outlierIndexes);
%---------------------------------------------------------------
% Fancy plots in the following section.  Delete if you don't need it.
% Show the original data and the absolute deviation
% subplot(2, 1, 1);
% bar(data);
% hold on;
% line(xlim, [meanValue, meanValue], 'Color', 'r', 'LineWidth', 2);
% grid on;
% title('Original Data', 'FontSize', fontSize);
% message = sprintf('Mean Value = %.2f', meanValue);
% text(3, 150, message, 'FontSize', 18, 'Color', 'r');
% subplot(2, 1, 2);
% bar(absoluteDeviation);
% grid on;
% title('Absolute Deviations', 'FontSize', fontSize);
% % Put a line for the mad.
% line(xlim, [mad, mad], 'Color', 'r', 'LineWidth', 2);
% message = sprintf('MAD Value = %.2f', mad);
% text(3, 50, message, 'FontSize', 18, 'Color', 'r');
% % Put a line for the mad.
% line(xlim, [thresholdValue, thresholdValue], 'Color', 'r', 'LineWidth', 2);
% message = sprintf('Outlier Threshold Value = %.2f', thresholdValue);
% text(3, 200, message, 'FontSize', 18, 'Color', 'r');
% % Set up figure properties:
% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% % Give a name to the title bar.
% set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off');

end