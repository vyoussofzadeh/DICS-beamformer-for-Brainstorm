
if exist(fullfile(outd.sub,['ic_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['ic_',subj,'.mat']));
else
    
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath = fullfile(outd.sub,['f_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath)
        else
            event = ft_read_event(datafile);
            clear val
            for i=2:length(event)
                val(i) = event(i).value;
            end
            val1 = unique(val);
            Evnt_IDs = mode(val);
            %         switch task
            %             case 1
            %                 Evnt_IDs =  min(val1);
            %         end
            disp('preprocessing ...');
            cfg = [];
            cfg.eventid = Evnt_IDs;
            cfg.epochtype = event(20).type;
            cfg.datafile  = datafile;
            cfg.hpfreq = 0.1;
            cfg.lpfreq = 40;
            [f_data] = vy_preprocess(cfg);
            [f_data] = vy_preprocess(cfg);
            disp('filtering was completed');
            save(savepath, 'f_data', '-v7.3');
        end
    end
    
    switch task
        case 2
            %% Adjust Offset Of Trial
            sampling_frequency = f_data.hdr.Fs;
            cfg = [];
            cfg.offset = -36 * sampling_frequency / 1000; % samples
            f_data = ft_redefinetrial(cfg, f_data);
    end
    
    %%
    %     cfg = [];
    %     cfg.resamplefs  = 250;
    %     cfg.demean      = 'no';
    %     cfg.detrend     = 'no';
    %     f_data = ft_resampledata(cfg, f_data);
    
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1
        
        savepath = fullfile(outd.sub,['a_',subj,'.mat']);
        cfg = [];
        cfg.pflag    = 2; % yes:1, No:2
        cfg.saveflag = 1; % yes:1, No:2
        cfg.savepath = savepath;
        cfg.latency  = [f_data.time{1}(1),f_data.time{1}(end)];%[-200,900]./1000;
        cfg.rejectpercentage = .95;
        [r_data,report] = vy_artifactreject(cfg, f_data);
        % disp('Bad data rejection was completed');
    end
    %% Bad trials & channels (Manuual)
    % clear r_data
    % cfg = [];
    % cfg.metric = 'zvalue';  % use by default zvalue method
    % cfg.latency = [-200,900];
    % cfg.layout   = lay;   % this allows for plotting individual trials
    % r_data   = ft_rejectvisual(cfg, f_data);
    
    %% Inspecting bad data
    % cfg = [];
    % cfg.viewmode = 'vertical';
    % cfg.continuous = 'no';
    % cfg.trials     = report.btrl;
    % cfg.channel   = report.bchan;
    % ft_databrowser(cfg,f_data);
    
    %% ICA cleaning
    warning off
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
        cfg.savefile = fullfile(cfg.savepath,['ic_',cfg.subj,'.mat']);
%         cln_data = vy_ica_cleaning_light(cfg, r_data);
        cln_data = vy_ica_cleaning_light2(cfg, r_data);
        disp('ICA cleaning was completed');
    end
end