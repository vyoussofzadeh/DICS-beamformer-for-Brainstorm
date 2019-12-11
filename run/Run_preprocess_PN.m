if exist(fullfile(outd.sub,['Scrambled_cln_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['Scrambled_cln_',subj,'.mat']));
    load(fullfile(outd.sub,['Images_cln_',subj,'.mat']));
else
    
    %     event = ft_read_event(datafile);
    %     clear val
    %     for i=1:length(event)
    %         val(i) = event(i).value;
    %     end
    %% Filteting, Event reading, Artifact rejecrtion
    savepath = fullfile(outd.sub,['full_f_',subj,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        disp('preprocessing ...');
        cfg = [];
        cfg.eventid = Evnt_IDs;
        cfg.epochtype = 'STI101';
        cfg.datafile  = datafile;
        cfg.hpfreq = 0.1;
        cfg.lpfreq = 40;
        [f_data, ecg_data] = vy_preprocess(cfg);
        disp('filtering was completed');
%         save(savepath, 'f_data', '-v7.3');
    end
    
    %% Bad trials & channels
    savepath = fullfile(outd.sub,['full_a_',subj,'.mat']);
    cfg = [];
    cfg.pflag = 2; % yes:1, No:2
    cfg.saveflag = 2; % yes:1, No:2
    cfg.savepath = savepath;
    cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];%[-200,900]./1000;
    [r_data,report] = vy_artifactreject(cfg, f_data);
    
    %% ICA cleaning
    cfg = [];
    cfg.savepath = fullfile(outd.sub,['full_ic_',subj,'.mat']);
    cfg.saveflag = 1;
    cfg.lay = lay;
    cfg.n   = 20;
    cln_data = vy_ica_cleaning_light(cfg, r_data);
    disp('ICA cleaning was completed');
end

%% SPLIT THE DATA
% 2: Scrambled images
savepath = fullfile(outd.sub,['Scrambled_cln_',subj,'.mat']);
if exist(savepath, 'file') == 2
    load(savepath)
else
    cfg = [];
    cfg.trials = find(cln_data.trialinfo == Evnt_IDs(1));
    Scrambled_cln_data = ft_selectdata(cfg, cln_data);
    save(savepath, 'Scrambled_cln_data', '-v7.3');
end

% 3: Images
savepath = fullfile(outd.sub,['Images_cln_',subj,'.mat']);
if exist(savepath, 'file') == 2
    load(savepath)
else
    cfg = [];
    cfg.trials = find(cln_data.trialinfo == Evnt_IDs(2));
    Images_cln_data = ft_selectdata(cfg, cln_data);
    save(savepath, 'Images_cln_data', '-v7.3');
end

