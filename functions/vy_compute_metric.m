function [level, info] = vy_compute_metric(cfg, data)
% SUBFUNCTION for ft_rejectvisual

% determine the initial selection of trials
ntrl = length(data.trial);
if isequal(cfg.trials, 'all') % support specification like 'all'
  cfg.trials = 1:ntrl;
end
trlsel = false(1, ntrl);
trlsel(cfg.trials) = true;

% determine the initial selection of channels
nchan = length(data.label);
cfg.channel = ft_channelselection(cfg.channel, data.label); % support specification like 'all'
chansel = false(1, nchan);
chansel(match_str(data.label, cfg.channel)) = true;

% compute the sampling frequency from the first two timepoints
fsample = 1/mean(diff(data.time{1}));

% select the specified latency window from the data
% here it is done BEFORE filtering and metric computation
for i=1:ntrl
  begsample = nearest(data.time{i}, cfg.latency(1));
  endsample = nearest(data.time{i}, cfg.latency(2));
  data.time{i} = data.time{i}(begsample:endsample);
  data.trial{i} = data.trial{i}(:, begsample:endsample);
end

% compute the offset from the time axes
offset = zeros(ntrl, 1);
for i=1:ntrl
  offset(i) = time2offset(data.time{i}, fsample);
end

% set up guidata info
info                = [];
info.data           = data;
info.cfg            = cfg;
info.metric         = cfg.metric;
info.previousmetric = 'none';
info.level          = nan(nchan, ntrl);
info.ntrl           = ntrl;
info.nchan          = nchan;
info.trlsel         = trlsel;
info.chansel        = chansel;
info.fsample        = fsample;
info.offset         = offset;
info.quit           = 0;


% update_log(cfg.output_box, 'Computing metric...');
% ft_progress('init', cfg.cfg.feedback, 'computing metric');
level = zeros(info.nchan, info.ntrl);
if strcmp(info.metric, 'zvalue') || strcmp(info.metric, 'maxzvalue')
  % cellmean and cellstd (see ft_denoise_pca) would work instead of for-loops, but they are too memory-intensive
  runsum = zeros(info.nchan, 1);
  runss  = zeros(info.nchan, 1);
  runnum = 0;
  for i=1:info.ntrl
    dat = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)),[]); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
    runsum = runsum + nansum(dat, 2);
    runss  = runss  + nansum(dat.^2, 2);
    runnum = runnum + sum(isfinite(dat), 2);
  end
  mval = runsum./runnum;
  sd   = sqrt(runss./runnum - (runsum./runnum).^2);
end
for i=1:info.ntrl
  ft_progress(i/info.ntrl, 'computing metric %d of %d\n', i, info.ntrl);
  dat = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)),[]); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
  switch info.metric
    case 'var'
      level(:, i) = nanstd(dat, [], 2).^2;
    case 'min'
      level(:, i) = nanmin(dat, [], 2);
    case 'max'
      level(:, i) = nanmax(dat, [], 2);
    case 'maxabs'
      level(:, i) = nanmax(abs(dat), [], 2);
    case 'range'
      level(:, i) = nanmax(dat, [], 2) - nanmin(dat, [], 2);
    case 'kurtosis'
      level(:, i) = kurtosis(dat, [], 2);
    case '1/var'
      level(:, i) = 1./(nanstd(dat, [], 2).^2);
    case 'zvalue'
      level(:, i) = nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2);
    case 'maxzvalue'
      level(:, i) = nanmax( ( dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , [], 2);
    otherwise
      ft_error('unsupported method');
  end
end