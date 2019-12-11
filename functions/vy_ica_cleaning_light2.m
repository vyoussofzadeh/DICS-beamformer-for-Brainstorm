function cln_data = vy_ica_cleaning_light2(cfg_main, data_in)

% satis = 0;
disp('ica cleaning ...');
if exist(cfg_main.savefile, 'file') && (cfg_main.overwrite) == 2
    load(cfg_main.savefile)
else
    
    savepath = fullfile(cfg_main.savepath,'ica');
    if exist(savepath, 'file') == 0, mkdir(savepath), end
    
    n = cfg_main.n; % ICs
    
    cfg = [];
    cfg.lay = cfg_main.lay;
    cfg.subj = cfg_main.subj;
    cfg.n = n;
    cfg.savefig = 1;
    cfg.allpath = cfg_main.allpath;
    comp = vy_ica(cfg, data_in);
    %     title(savepath)
    cfg = [];
    cfg.updatesens = 'no';
    
    %%
    cfg = [];
    cfg.pflag    = 2; % yes:1, No:2
    cfg.saveflag = 2; % yes:1, No:2
    cfg.savepath = [];
    cfg.latency  = [comp.time{1}(1),comp.time{1}(end)];%[-200,900]./1000;
    cfg.rejectpercentage = .95;
    [IC_data,report] = vy_artifactreject(cfg, comp);
    save('ica/ICAreport.mat','report');
    disp(report)
    
    %% Rejecting bad trials, identified by the ICA
    trials = find(~ismember(1:length(comp.trial),report.btrl));
    cfg = [];
    cfg.trials = trials;
    data_in = ft_selectdata(cfg, data_in);
    
    %%
    %     cfg = [];
    %     cfg.viewmode = 'component';
    %     cfg.layout = cfg_main.lay;
    %     ft_databrowser(cfg, r_data);
    %     colormap(brewermap(256, '*RdYlBu'));
    %     set(gcf, 'Position', [600   600   700   500]);
    
    %% Rejecting bad ICAs and prjecting back into data space
    disp('=============================')
    disp('suggested')
    disp(report.bchan)
    cfg = [];
    if cfg_main.select == 1
        bic = input(['Select bad ICs for ' cfg_main.subj,':']);
        cfg.component = comp.label(bic);
    else
        cfg.component =  report.bchan;
    end
    cfg.updatesens = 'no';
    cfg.trials = trials;
    cln_data = ft_rejectcomponent(cfg, comp, data_in);
    close all
    if cfg_main.saveflag == 1
        save(cfg_main.savefile, 'cln_data', '-v7.3');
        textfile_rej = 'ica/selected_badICs';
        badICs = cell2table(cfg.component);
        badICs.Properties.VariableNames{'Var1'} = 'bICAs';
        writetable(badICs,textfile_rej,'Delimiter',' ');        
    end
    
    %% back to new ft
    restoredefaultpath
    addpath((cfg_main.allpath.ft_path));
    ft_defaults
    addpath(genpath(cfg_main.allpath.hcp_path));
    addpath(genpath(cfg_main.allpath.cd_org));
    addpath(genpath(cfg_main.allpath.exfig_path));
    
end

