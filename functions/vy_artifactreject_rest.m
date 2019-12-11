function [r_data,report] = vy_artifactreject_rest(cfg_main, dat)

disp('Bad channels, bad trials ...');
if exist(cfg_main.savepath, 'file') == 2
    load(cfg_main.savepath)
else
    
    %% kurtosis
    cfg = [];
    cfg.trials = 'all';
    cfg.metric = 'kurtosis';
    cfg.channel = 'all';
    cfg.latency = cfg_main.latency;
    [level,info] = vy_compute_metric(cfg,dat); metric.kurt = level;
    info.pflag = cfg_main.pflag;
    [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level);
    
    thresh.(cfg.metric) = 0.95.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
    thresh.(cfg.metric) = 0.95.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel
    
    %% zvalue
    cfg = [];
    cfg.trials = 'all';
    cfg.metric = 'zvalue';
    cfg.channel = 'all';
    cfg.latency = cfg_main.latency;
    [level,info] = vy_compute_metric(cfg,dat); metric.kurt = level;
    info.pflag = cfg_main.pflag;
    [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level);
    
    thresh.(cfg.metric) = 0.95.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
    thresh.(cfg.metric) = 0.95.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel
    
    %% Var
    cfg = [];
    cfg.trials = 'all';
    cfg.metric = 'var';
    cfg.channel = 'all';
    cfg.latency = cfg_main.latency;
    [level,info] = vy_compute_metric(cfg,dat); metric.kurt = level;
    info.pflag = cfg_main.pflag;
    [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level);
    
    thresh.(cfg.metric) = 0.95.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
    thresh.(cfg.metric) = 0.95.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel
    
    %%
%     btrl_all = unique([btrl.kurtosis,btrl.var,btrl.zvalue]);
    bch_all = unique([bch.kurtosis;bch.var;bch.zvalue]);
%     disp('Bad trials:')
%     disp(btrl_all);
    disp('Bad channels:')
    for i=1:length(bch_all)
        bch_all_label_disp{i,:} = dat.label{bch_all(i)};
        bch_all_label{i,:} = ['-',dat.label{bch_all(i)}];
    end
    disp(bch_all_label_disp);
    
    %% Removing bad channels/trials
%     cfg = [];
%     cfg.trials = find(~ismember(1:length(dat.trial),btrl_all));
%     dat = ft_selectdata(cfg, dat);
    
    if length(bch_all) < 10
        cfg = [];
        cfg.channel = ['all';bch_all_label];
        dat = ft_selectdata(cfg, dat);
    else
        warning('bad chans were more than 8 (arbitary), have a look at the data');
    end
    
%     report.btrl = btrl_all; 
    report.bchan = bch_all_label_disp;   
    r_data = dat;
    if cfg_main.saveflag ==1
        save(cfg_main.savepath, 'r_data', 'report','-v7.3');
    end
end




%%
% cfg = [];
% cfg.trials = 'all';
% cfg.metric = 'kurtosis';
% cfg.channel = 'all';
% cfg.latency = [-400,900];
% [level,info] = vy_compute_metric(cfg,f_data2); metric.kurt = level;
% pflag = 1;
% [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level,pflag);
