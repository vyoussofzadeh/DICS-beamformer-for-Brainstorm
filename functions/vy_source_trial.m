function [s_data, s_data2] = vy_source_trial(data, grid, vol)

% lambda = 5; % percent
% data_cov = squeeze(mean(data.app.cov,1));
% s = svd(data_cov);
% figure,plot(s,'.');
% hold on
% reg = lambda/100*mean(s);
% plot([0,length(s)],[reg,reg],'r')

%
% lambda = .1;

% lambda_bsl = 5; % percent
% data_cov = squeeze(mean(data.bsl.cov,1));
% s = svd(data_cov);
% figure,
% plot(s,'.r');
% hold on
% reg = lambda_bsl/100*mean(s);
% plot([0,length(s)],[reg,reg],'r')
% 
% 
% lambda_pst = 5; % percent
% data_cov = squeeze(mean(data.pst.cov,1));
% s = svd(data_cov);
% % figure,
% plot(s,'.b');
% % hold on
% reg = lambda_pst/100*mean(s);
% plot([0,length(s)],[reg,reg],'b')

%%

% create spatial filter using the lcmv beamformer
cfg                  = [];
cfg.method           = 'lcmv';
cfg.grid             = grid; % leadfield, which has the grid information
cfg.headmodel         = vol; % volume conduction model (headmodel)
cfg.keepfilter       = 'yes';
cfg.lcmv.keepfilter  = 'yes';
cfg.keeptrials       = 'yes';
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.lambda      = '5%';

sourceAll = ft_sourceanalysis(cfg, data.app);
cfg.grid.filter = sourceAll.avg.filter;
s_data.bsl = ft_sourceanalysis(cfg, data.bsl);
s_data.pst = ft_sourceanalysis(cfg, data.pst);

cfg.rawtrial = 'yes';
s_data2.bsl  = ft_sourcedescriptives([], s_data.bsl); % to get the neural-activity-index
s_data2.pst  = ft_sourcedescriptives([], s_data.pst); % to get the neural-activity-index
