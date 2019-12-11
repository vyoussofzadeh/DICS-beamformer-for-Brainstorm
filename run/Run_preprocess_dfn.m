if exist(fullfile(outd.sub,['full_ic_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['full_ic_',subj,'.mat'])); cln.left = cln_data;
    load(fullfile(outd.sub,['full_ic_',subj,'.mat'])); cln.right = cln_data;
else
    
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath = fullfile(outd.sub,['full_f_',subj,'.mat']);
%         savepath2 = fullfile(outd.sub,['f_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath)
%             load(savepath2)
        else
            cfg = [];
            cfg.eventid = [1,2];
            cfg.epochtype = epoch_type;
            cfg.datafile  = datafile;
            cfg.hpfreq = 0.1;
            cfg.lpfreq = 40;
            f_data = vy_preprocess_dfn(cfg);
%             cfg.eventid = 2;
%             [f_nosie_data] = vy_preprocess_motor(cfg);
            disp('filtering was completed');
            save(savepath, 'f_data', '-v7.3');
%             save(savepath2, 'fr_data', '-v7.3');
        end
    end
    
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1       
%         savepath1 = fullfile(outd.sub,['al_',subj,'.mat']);
        savepath = fullfile(outd.sub,['full_a_',subj,'.mat']);
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 2; % yes:1, No:2
        cfg.savepath = savepath;
        cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)]; %[-200,900]./1000;
        cfg.rejectpercentage = .95;
        [r_data,report] = vy_artifactreject(cfg, f_data);
%         cfg.latency = [fr_data.time{1}(1),fr_data.time{1}(end)]; %[-200,900]./1000;
%         cfg.savepath = savepath2;
%         [rr_data,report_r] = vy_artifactreject(cfg, fr_data);
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
        cfg.savefile = fullfile(cfg.savepath,['full_ic_',cfg.subj,'.mat']);
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        cln_data = vy_ica_cleaning_light2(cfg, r_data);
%         cfg.savefile = fullfile(cfg.savepath,['icr_',cfg.subj,'.mat']);
%         clnr_data = vy_ica_cleaning_light2(cfg, rr_data);
        disp('ICA cleaning was completed');
    end
end

cfg = [];
cfg.trials = find(cln_data.trialinfo == 1);
data_task = ft_selectdata(cfg, cln_data);

cfg.trials = find(cln_data.trialinfo == 2);
data_noise = ft_selectdata(cfg, cln_data);

% data_all = ft_appenddata([], data_noise,data_task);
% save('ic_carver_brian.mat', 'cln_data', '-v7.3');