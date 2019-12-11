
if exist(fullfile(outd.sub,['icl_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['icl_',subj,'.mat'])); cln.left = cln_data;
    load(fullfile(outd.sub,['icr_',subj,'.mat'])); cln.right = cln_data;
else
    
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath1 = fullfile(outd.sub,['fl_',subj,'.mat']);
        savepath2 = fullfile(outd.sub,['fr_',subj,'.mat']);
        if exist(savepath1, 'file') == 2
            load(savepath1)
            load(savepath2)
        else
            cfg = [];
            cfg.eventid = 16;
            cfg.epochtype = epoch_type;
            cfg.datafile  = datafile;
            cfg.hpfreq = 0.1;
            cfg.lpfreq = 40;
            [fl_data] = vy_preprocess_motor(cfg);
            cfg.eventid = 32;
            [fr_data] = vy_preprocess_motor(cfg);
            disp('filtering was completed');
            save(savepath1, 'fl_data', '-v7.3');
            save(savepath2, 'fr_data', '-v7.3');
        end
    end
    
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1       
        savepath1 = fullfile(outd.sub,['al_',subj,'.mat']);
        savepath2 = fullfile(outd.sub,['ar_',subj,'.mat']);
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 2; % yes:1, No:2
        cfg.savepath = savepath1;
        cfg.latency = [fl_data.time{1}(1),fl_data.time{1}(end)]; %[-200,900]./1000;
        cfg.rejectpercentage = .95;
        [rl_data,report_l] = vy_artifactreject(cfg, fl_data);
        cfg.latency = [fr_data.time{1}(1),fr_data.time{1}(end)]; %[-200,900]./1000;
        cfg.savepath = savepath2;
        [rr_data,report_r] = vy_artifactreject(cfg, fr_data);
        disp('Bad data rejection was completed');
    end
    
    %% ICA cleaning
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        %         cfg.savepath = fullfile(outd.sub,['ic_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        cfg.savefile = fullfile(cfg.savepath,['icl_',cfg.subj,'.mat']);
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        clnl_data = vy_ica_cleaning_light2(cfg, rl_data);
        cfg.savefile = fullfile(cfg.savepath,['icr_',cfg.subj,'.mat']);
        clnr_data = vy_ica_cleaning_light2(cfg, rr_data);
        disp('ICA cleaning was completed');
    end
end