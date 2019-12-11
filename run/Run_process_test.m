clear; clc, close('all'); warning off

%% Initial settings
set(0,'DefaultFigureWindowStyle','docked')

cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%-- Input dir
indir = '/data/MEG/Vahab/test_data';
%-- Output dir
outdir = '/data/MEG/Vahab/test_data/processed';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.anatomy = 1;     % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis

%%
disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');

switch task
    case 1
        % - Auditory definition naming
        tag = 'DFN'; tag1 = 'dfn';
        Evnt_IDs = 1; % questions
    case 2
        % - Visual picture naming
        tag = 'PN';
        Evnt_IDs = 3; % 3: images, 2: scrambled images
end

%% analysis flag
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis

%%
disp('1: 2019')
disp('2: 2018');
disp('3: 2017');
disp('4: 2016')
disp('5: 2015');
disp('6: 2014');
disp('7: 2013')
disp('8: 2012');
disp('9: 2011');
disp('10: Older');
year = input('Year data were acquired: ');

clear ytag;
switch year
    case 1
        ytag = {'19'};
    case 2
        ytag = {'18'};
    case 3
        ytag = {'17'};
    case 4
        ytag = {'16'};
    case 5
        ytag = {'up'};
    case 6
        ytag = {'up'};
    case 7
        ytag = {'13'};
    case 8
        ytag = {'12'};
    case 9
        ytag = {'11'};
end
disp('============');

%% All data
% clear datafolder datafile
% datafile1 = [];
% d = rdir([indir,['/**/','sss','/*',tag,'*/*raw_tsss.fif']]);
% d = rdir([indir,['/**/','/*',tag,'*/*raw_tsss.fif']]);
%
% for i=1:length(d)
%     [pathstr, name] = fileparts(d(i).name);
%     datafolder{i} = pathstr;
%     datafile{i} = d(i).name;
% end
% datafile1 = vertcat(datafile1,datafile);
% datafile1 = datafile1';
% disp(datafile1)

%% Listing data
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
% Per year
clear datafolder datafile subj_all
datafile1 = [];
d = rdir([indir,['/**/','/*',tag,'*/*raw_tsss.fif']]);
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj_all{i} = [num2str(i), ': ', datafile{i}(Index(5)+1:Index(6)-1)];
end
datafile1 = vertcat(datafile1,datafile);
datafile1 = datafile1';
disp(datafile1)
disp('============');

%%
disp('1: choose specific subject');
disp('2: do all');
subsel = input('?');
switch subsel
    case 1
        disp('Subjects')
        disp(subj_all')
        subsel1 = input('enter subject number?');
end
disp('============');

%%
epoch_type = 'STI101';

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);
disp('============');

%%
clear datafile2
switch subsel
    case 1
        datafile2{1} = datafile1{subsel1}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    case 2
        datafile2 = datafile1;
end

%%
for i = 1:size(datafile2,1)
    
    datafile = datafile2{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(5)+1:Index(6)-1);
    Date  = datafile(Index(6)+1:Index(7)-1);
    disp(datafile)
    disp(['subj:',subj])
    disp(['Date:',Date])
    disp('============');
    
    %-elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    disp('============');
    
    %%
    if year>=4
        yttag = 'older';
    else
        yttag = ytag{1};
    end
    outd.sub = fullfile(outdir,'ft_process',yttag, subj, tag);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    disp(['outputdir:',outd.sub])
    disp('============');
    
    %% Preprocesssing
    ic_selection = 1; % 1: manual, 2: automated
    Run_preprocess   
    
    %% Notch filtering of 30Hz
    if flag.notch ==1
        Run_notch
    end
    
    %% FFT & TFR
    datain = cln_data;
    if flag.freq == 1
        stag = 'tsk_baseline';
        Run_freq
        disp(['time_of_interest:',num2str(time_of_interest),'sec']);
        disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
        L = 0.3;
    end
    
    %%
    %     cfg = [];
    %     cfg.savefile = [];
    %     cfg.saveflag = 2;
    %     cfg.foilim = [2 40];
    %     cfg.plotflag  = 1;
    %     cfg.tapsmofrq = 5;
    %     cfg.taper     = 'hanning';
    %     vy_fft(cfg, cln_data);
    %     grid on
    %     grid minor
    %     title('Before band-stop filtering');
    %
    %     %%
    %     fsb = input('Enter the sop-band frequency?');
    %     cfg = [];
    %     cfg.bsfilter = 'yes';
    %     %     cfg.bsfreq = [29 32]; % or whatever you deem appropriate
    %     cfg.bsfreq = [fsb-1 fsb+1]; % or whatever you deem appropriate
    %     %     cfg.bsfreq = [8 12;29 32]; % or whatever you deem appropriate
    %     cln_data = ft_preprocessing(cfg, cln_data);
    %     %     cfg.bsfreq = [2 12]; % or whatever you deem appropriate
    %
    %     cfg = [];
    %     cfg.savefile = [];
    %     cfg.saveflag = 2;
    %     cfg.foilim = [2 40];
    %     cfg.plotflag  = 1;
    %     cfg.tapsmofrq = 5;
    %     cfg.taper     = 'hanning';
    %     vy_fft(cfg, cln_data);
    %     grid on
    %     grid minor
    %     title('After band-stop filtering');
    
    %% FFT & TFR
    if flag.freq == 1
        stag = 'tsk_baseline'; datain = cln_data;
        Run_freq
        disp(['time_of_interest:',num2str(time_of_interest),'sec']);
        disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
        L = 0.3;
    end
    
    %% Timelock
    if flag.time == 1
        %--setting baseline interval
        toi(1,:) = [-0.3,0];
        
        %-- setting the post-stim interval
        disp(['suggested time interval:',num2str(time_of_interest), '+-', num2str(L),' Sec']);
        %     disp('Yes: 1, No: 2');
        %     tfa = input('Is it OK to proceed?');
        %     disp('the following time was selected');
        %        if (time_of_interest-L) > 0.4 && (time_of_interest+L) < 2.3
        if (time_of_interest-L) >= 0.4 && (time_of_interest+L) < 2.3
            toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
        else
            switch task
                case 1
                    %                         toi = [-0.3,0;1.1,1.7]; % Best of DN
                    toi(2,:) = [0.8,1.5]; % Best of DN
                case 2
                    toi(2,:) = [0.4,1.2]; % Best for PN, left IFG
                    toi(2,:) = [0.4,1.6]; % Best for PN, left IFG
                    toi(2,:) = [0.7,1.6]; % Best for PN, left IFG
                    toi(2,:) = [1,1.6]; % Best for PN, left IFG
            end
        end
        disp(['[',num2str(toi(1,1)),',',num2str(toi(1,2)),'] sec interval was selected as bsl']);
        disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected as pst']);
        Run_time
    end
%     disp(['suggested time interval:',num2str(time_of_interest), '+-', num2str(L),' Sec']);
%     disp('Yes: 1, No: 2');
%     tfa = input('Is it OK to proceed?');
%     disp('the following time was selected');
%     
%     if flag.time == 1
%         switch tfa
%             case 1
%                 toi(1,:) = [-0.3,0];
%                 toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
%             case 2
%                 switch task
%                     case 1
%                         toi = [-0.3,0;1.1,1.7]; % Best of DN
%                         %                         toi = [-0.3,0;.8,1.5]; % Best of DN
%                     case 2
%                         toi = [-0.3,0;0.4,1.2]; % Best for PN, left IFG
%                         %                         toi = [-0.3,0;0.8,1.5]; % Best of DN
%                 end
%         end
%         disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected']);
%         pause(2)
%         Run_time
%     end
    
    %% Grand Mean
    if flag.gave == 1
        Run_grandmean
    end
    %% Source analysis
    outputmridir = fullfile(outdir,'ft_process',yttag, subj,'anat'); % output dir
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    
    %% Processing
    disp('1: Surface-based')
    disp('2: Volumetric');
    analysis = input('Eneter the analysis: ');
%     analysis = 2;
    
    switch analysis
        case 1
            mtag = 'source_surf';
            disp('1: SPM source analysis (surface + BF)');
            disp('2: Export ft to Brainstorm');
            disp('3: Source-based bf');
            method = input('Method: ');
        case 2
            % end
            disp('1: LCMV')
            disp('2: Network+Connectvity');
            disp('3: DICS Source');
            method = input('Method: ');
%             method = 3;
            %-
            %             disp('1: Low-res grid')
            %             disp('2: High-res grid')
            %             meshgrid = input('Mesh grid: ');
            meshgrid = 1;
    end
    
    switch analysis
        case 1
            % Surface-based analysis
            Run_surfacebased
        case 2
            % Volumetric-based analysis
            anatomy_check_flag = 2;
            Run_volumetric
    end
    disp('============');
    
    %%
    pause
    close all
    disp([datafile,' ,was completed'])
    
end




