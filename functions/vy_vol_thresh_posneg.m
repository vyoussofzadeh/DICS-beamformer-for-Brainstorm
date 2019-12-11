function s_vol = vy_vol_thresh_posneg(s_vol,projthresh,posneg)

% set threshold
% projthresh = 0.5;
val = s_vol.anatomy;
% val(abs(val) < projthresh*max(abs(val(:)))) = NaN;
if posneg == 1
    val(val < projthresh) = NaN;
elseif posneg == -1
    val(val > projthresh) = NaN; % keep negative effects
end
val = val./max(abs(val(:)));
s_vol.anatomy = val;
