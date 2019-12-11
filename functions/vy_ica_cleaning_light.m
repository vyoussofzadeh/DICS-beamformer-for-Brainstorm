function cln_data = vy_ica_cleaning_light(cfg_main, f_data)

n = cfg_main.n; % ICs

% satis = 0;
disp('ica cleaning ...');
if exist(cfg_main.savepath, 'file') == 2
    load(cfg_main.savepath)
else
    
    cfg = [];
    cfg.lay = cfg_main.lay;
    cfg.subj = cfg_main.subj;
    cfg.n = n;
    cfg.allpath = cfg_main.allpath;
    comp = vy_ica(cfg, f_data);
    title(savepath)
    cfg = [];
    cfg.updatesens = 'no';
    disp('=============================')
    bic = input(['Select bad ICs for slected data for ' cfg_main.subj,':']);
    cfg.component = comp.label(bic);
    cln_data = ft_rejectcomponent(cfg, comp, f_data);
    close all
    if cfg_main.saveflag == 1
        save(cfg_main.savepath, 'cln_data', '-v7.3');
    end
    
    %% back to new ft
    restoredefaultpath
    addpath((cfg_main.allpath.ft_path));
    ft_defaults
    addpath(genpath(cfg_main.allpath.hcp_path));
    addpath(genpath(cfg_main.allpath.cd_org));
    addpath(genpath(cfg_main.allpath.exfig_path));

end

