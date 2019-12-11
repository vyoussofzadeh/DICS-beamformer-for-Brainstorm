function [s_data, s_data2] = vy_source_freq(cfg_main, data)

mtag = cfg_main.mtag;

switch mtag
    
    case {'dics','dics_ratio'}
        %         cfg              = [];
        %         cfg.method       = 'dics';
        %         cfg.frequency    = data.app.freq;
        %         cfg.grid         = grid;
        %         cfg.headmodel    = vol;
        %         cfg.dics.lambda       = '5%';
        %         cfg.dics.keepfilter   = 'yes';
        %         %         cfg.dics.realfilter   = 'yes';
        %         cfg.dics.fixedori     = 'yes';
        %         sourceAll = ft_sourceanalysis(cfg, data.app);
        %         cfg.grid.filter = sourceAll.avg.filter;
        %         cfg.dics.lambda       = '5%';
        %         s_data.bsl = ft_sourceanalysis(cfg, data.bsl);
        %         s_data.pst = ft_sourceanalysis(cfg, data.pst);
        
        cfg = [];
        cfg.method = 'dics';
        cfg.dics.lambda = '100%';
        cfg.frequency    = data.app.freq;
        cfg.grid = cfg_main.grid;
        cfg.headmodel = cfg_main.headmodel;
        cfg.dics.keepfilter = 'yes';
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        sourceavg = ft_sourceanalysis(cfg, data.app);
        
        cfg = [];
        cfg.method = 'dics';
        %         cfg.dics.lambda = '5%';
        cfg.grid = cfg_main.grid;
        cfg.grid.filter = sourceavg.avg.filter;
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.headmodel = cfg_main.headmodel;
        s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst      = ft_sourceanalysis(cfg, data.pst);
        
        %         cfg = [];
        %         cfg.method = 'dics';
        %         cfg.dics.lambda = '5%';
        %         cfg.frequency    = data.app.freq;
        %         cfg.grid = cfg_main.grid;
        %         cfg.headmodel = cfg_main.headmodel;
        %         cfg.dics.keepfilter = 'yes';
        %         cfg.dics.realfilter   = 'yes';
        %         cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        %         sourceavg = ft_sourceanalysis(cfg, data.app);
        %
        %         cfg = [];
        %         cfg.method = 'dics';
        %         cfg.dics.lambda = '5%';
        %         cfg.grid = cfg_main.grid;
        %         cfg.grid.filter = sourceavg.avg.filter;
        %         cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        %         cfg.headmodel = cfg_main.headmodel;
        %         s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
        %         s_data.pst      = ft_sourceanalysis(cfg, data.pst);
        %
        
    case 'dics_fs'
        
        %         cfg = [];
        %         cfg.method = 'dics';
        %         cfg.dics.lambda = '100%';
        %         %         cfg.frequency    = data.app.freq;
        %         cfg.sourcemodel              = cfg_main.sourcemodel;
        %         %         cfg.grid.leadfield    = cfg_main.leadfield.leadfield;
        %         cfg.headmodel         = cfg_main.headmodel;
        %         cfg.dics.keepfilter   = 'yes';
        %         cfg.dics.fixedori     = 'yes'; % project on axis of most variance using SVD
        %         sourceavg = ft_sourceanalysis(cfg, data.app);
        %         %         sourceavg = ft_sourceanalysis(cfg, ft_checkdata(data.app,'cmbrepresentation','fullfast')); % trick to speed up the computation
        %
        %         cfg = [];
        %         cfg.method = 'dics';
        %         %         cfg.dics.lambda = '0%';
        %         %         cfg.leadfield             = cfg_main.leadfield;
        %         cfg.sourcemodel              = cfg_main.sourcemodel;
        %         %         cfg.grid              = cfg_main.sourcemodel;
        %         %         cfg.grid.leadfield    = cfg_main.leadfield.leadfield;
        %         cfg.headmodel        = cfg_main.headmodel;
        %         cfg.grid.filter = sourceavg.avg.filter;
        %         cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        %         s_data.bsl           = ft_sourceanalysis(cfg, data.bsl);
        %         s_data.pst           = ft_sourceanalysis(cfg, data.pst);
        %
        %         cfg = [];
        %         cfg.method = 'dics';
        %         cfg.dics.lambda = '100%';
        % %         cfg.frequency    = data.app.freq;
        %         cfg.sourcemodel              = cfg_main.sourcemodel;
        % %         cfg.grid.leadfield    = cfg_main.leadfield.leadfield;
        %         cfg.headmodel         = cfg_main.headmodel;
        %         cfg.dics.keepfilter   = 'yes';
        %         cfg.dics.fixedori     = 'yes'; % project on axis of most variance using SVD
        %         sourceavg = ft_sourceanalysis(cfg, data.app);
        % %         sourceavg = ft_sourceanalysis(cfg, ft_checkdata(data.app,'cmbrepresentation','fullfast')); % trick to speed up the computation
        %
        %         cfg = [];
        %         cfg.method = 'dics';
        % %         cfg.dics.lambda = '0%';
        %         cfg.grid             = cfg_main.grid;
        % %         cfg.sourcemodel              = cfg_main.grid;
        % %         cfg.grid              = cfg_main.sourcemodel;
        % %         cfg.grid.leadfield    = cfg_main.leadfield.leadfield;
        %         cfg.headmodel        = cfg_main.headmodel;
        %         cfg.grid.filter = sourceavg.avg.filter;
        %         cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        %         s_data.bsl           = ft_sourceanalysis(cfg, data.bsl);
        %         s_data.pst           = ft_sourceanalysis(cfg, data.pst);
        
        
        %         cfg = [];
        %         cfg.method = 'dics';
        %         cfg.dics.lambda = '10%';
        %         %         cfg.frequency    = data.app.freq;
        %         cfg.grid              = cfg_main.grid;
        %         %         cfg.grid.leadfield    = cfg_main.leadfield.leadfield;
        %         cfg.headmodel         = cfg_main.headmodel;
        %         cfg.dics.keepfilter   = 'yes';
        %         cfg.dics.fixedori     = 'yes'; % project on axis of most variance using SVD
        %         sourceavg = ft_sourceanalysis(cfg, data.app);
        %         %         sourceavg = ft_sourceanalysis(cfg, ft_checkdata(data.app,'cmbrepresentation','fullfast')); % trick to speed up the computation
        %
        %         cfg = [];
        %         cfg.method = 'dics';
        %         cfg.dics.lambda = '5%';
        %         cfg.grid             = cfg_main.grid;
        %         %         cfg.sourcemodel              = cfg_main.grid;
        %         %         cfg.grid              = cfg_main.sourcemodel;
        %         %         cfg.grid.leadfield    = cfg_main.leadfield.leadfield;
        %         cfg.headmodel        = cfg_main.headmodel;
        %         cfg.grid.filter = sourceavg.avg.filter;
        %         cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        %         s_data.bsl           = ft_sourceanalysis(cfg, data.bsl);
        %         s_data.pst           = ft_sourceanalysis(cfg, data.pst);
        
        cfg = [];
        cfg.method = 'dics';
        cfg.dics.lambda = '100%';
        %         cfg.frequency    = data.app.freq;
        cfg.grid              = cfg_main.grid;
        %         cfg.grid.leadfield    = cfg_main.leadfield.leadfield;
        cfg.headmodel         = cfg_main.headmodel;
        cfg.dics.keepfilter = 'yes';
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        sourceavg = ft_sourceanalysis(cfg, data.app);
        %         sourceavg = ft_sourceanalysis(cfg, ft_checkdata(data.app,'cmbrepresentation','fullfast')); % trick to speed up the computation
        
        cfg = [];
        cfg.method = 'dics';
%         cfg.dics.lambda = '%';
        cfg.grid       = cfg_main.grid;
        %         cfg.grid.leadfield    = cfg_main.leadfield.leadfield;
        cfg.headmodel         = cfg_main.headmodel;
        cfg.grid.filter = sourceavg.avg.filter;
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst      = ft_sourceanalysis(cfg, data.pst);
        
        
    case 'dics_stat'
        
        cfg = [];
        cfg.method = 'dics';
        cfg.dics.lambda = '5%';
        cfg.frequency    = data.app.freq;
        % cfg.grid = grid;
        cfg.grid.pos = cfg_main.grid.pos(cfg_main.grid.inside,:);
        cfg.headmodel = cfg_main.headmodel;
        cfg.dics.keepfilter = 'yes';
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.channel = data.bsl.label;
        sourceavg = ft_sourceanalysis(cfg, data.app);
        
        cfg = [];
        cfg.method = 'dics';
        cfg.dics.lambda = '5%';
        cfg.grid.pos = cfg_main.grid.pos(cfg_main.grid.inside,:);
        cfg.grid.filter = sourceavg.avg.filter;
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.rawtrial = 'yes';
        cfg.headmodel = cfg_main.headmodel;
        s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst      = ft_sourceanalysis(cfg, data.pst);
        
    case 'pcc'
        
        cfg              = [];
        cfg.method       = 'pcc';
        cfg.frequency    = data.app.freq;
        cfg.grid         = grid;
        cfg.headmodel    = vol;
        cfg.pcc.lambda       = '5%';
        cfg.pcc.keepfilter   = 'yes';
        cfg.pcc.projectnoise  = 'yes';
        cfg.pcc.fixedori     = 'yes';
        sourceAll = ft_sourceanalysis(cfg, data.app);
        cfg.grid.filter = sourceAll.avg.filter;
        s_data.bsl = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst = ft_sourceanalysis(cfg, data.pst);
end

cfg = [];
% cfg.projectmom  = 'yes';
s_data2.bsl  = ft_sourcedescriptives(cfg, s_data.bsl); % to get the neural-activity-index
s_data2.pst  = ft_sourcedescriptives(cfg, s_data.pst); % to get the neural-activity-index




