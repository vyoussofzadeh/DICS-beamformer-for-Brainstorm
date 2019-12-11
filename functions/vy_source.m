function s_data = vy_source(cfg_main, data)

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
switch cfg_main.mtd
    
    case 'lcmv'
        % create spatial filter using the lcmv beamformer
        cfg                  = [];
        cfg.method           = 'lcmv';
        cfg.grid             = cfg_main.grid; % leadfield, which has the grid information
        cfg.headmodel         = cfg_main.headmodel; % volume conduction model (headmodel)
        cfg.keepfilter       = 'yes';
        cfg.lcmv.keepfilter  = 'yes';
        cfg.keeptrials       = 'yes';
        cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
        %         cfg.lcmv.lambda      = '5%';
        %         cfg.lcmv.lambda      = '100%';
        cfg.lcmv.lambda      = '0.1%';
        
        sourceAll = ft_sourceanalysis(cfg, data.app);
        
        % now applying the precomputed filters to pre and post intervals
        cfg                  = [];
        cfg.method           = 'lcmv';
        cfg.grid             = cfg_main.grid; % leadfield, which has the grid information
        cfg.headmodel        = cfg_main.headmodel; % volume conduction model (headmodel)
        cfg.grid.filter = sourceAll.avg.filter;
        s_data.bsl = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst = ft_sourceanalysis(cfg, data.pst);
        
        %         s_data2.bsl  = ft_sourcedescriptives([], s_data.bsl); % to get the neural-activity-index
        %         s_data2.pst  = ft_sourcedescriptives([], s_data.pst); % to get the neural-activity-index
        
    case 'lcmv_stat'
        
        cfg = [];
        cfg.method = 'lcmv';
        cfg.lcmv.lambda = '5%';
        % cfg.grid = grid;
        cfg.grid.pos = cfg_main.grid.pos(cfg_main.grid.inside,:);
        cfg.headmodel = cfg_main.headmodel;
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
        sourceavg = ft_sourceanalysis(cfg, data.app);
        
        cfg = [];
        cfg.method = 'lcmv';
        cfg.lcmv.lambda = '5%';
        cfg.grid.pos = cfg_main.grid.pos(cfg_main.grid.inside,:);
        cfg.grid.filter = sourceavg.avg.filter;
        cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.rawtrial = 'yes';
        cfg.headmodel = cfg_main.headmodel;
        s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst      = ft_sourceanalysis(cfg, data.pst);
        
    case 'sam'
        % create spatial filter using the lcmv beamformer
        cfg                  = [];
        cfg.method           = 'sam';
        cfg.grid             = cfg_main.grid; % leadfield, which has the grid information
        cfg.headmodel         = cfg_main.headmodel; % volume conduction model (headmodel)
        cfg.keepfilter       = 'yes';
        cfg.sam.keepfilter  = 'yes';
        cfg.keeptrials       = 'yes';
        %         cfg.sam.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.sam.lambda      = '5%';
        %         cfg.lcmv.lambda      = '100%';
        
        sourceAll = ft_sourceanalysis(cfg, data.app);
        
        % now applying the precomputed filters to pre and post intervals
        cfg                  = [];
        cfg.method           = 'lcmv';
        cfg.grid             = cfg_main.grid; % leadfield, which has the grid information
        cfg.headmodel        = cfg_main.headmodel; % volume conduction model (headmodel)
        cfg.grid.filter = sourceAll.avg.filter;
        s_data.bsl = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst = ft_sourceanalysis(cfg, data.pst);
        
        
end