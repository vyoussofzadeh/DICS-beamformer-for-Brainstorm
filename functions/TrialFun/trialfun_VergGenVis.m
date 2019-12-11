function    trl = trialfun_VergGenVis( cfg )

% This is the Trial Definition Function for Visual Object Naming [ Object(1.5s) + Cross(1.5s) + Blank (3.0s)].
% Subject is instructed to say object name aloud when see the Cross.
% The trigger is derived from the PhotoDiode signal of the trigger channle
%__________________________________________________________________________
%% Input:
% cfg: Structure with fields that contain parameters for extracting trials
% cfg.datafile: Char string represents filename of raw input MEG data
% cfg.trialdef.stimulCode: Respresting the desired stimulus type to be extracted from the data with following integer values:
% cfg.trialdef.stimulCode = [1,2,3];
% % trigger_name = {'Noune','+','Response'};
% % trigger_code = { 2176 , 2112 , 2176 };
% % 'All': Starts with Object and long enough until end of Naming

% cfg.trialdef.prestimTime: Desired time interval (in seconds) before Trigger onset. 
%       Must be positive. In 'blocks' mode, prestimTime w.r.t. first event.
% cfg.trialdef.poststimTime: Desired time interval (in seconds) after Trigger onset.
%       Must be positive. In 'blocks' mode, poststimTime w.r.t. last event.
% cfg.trialdef.plotresults = 'yes' or 'no' (default = 'no')
%__________________________________________________________________________
%% Output:
% trl: Ntrials-by-3 matrix with the trial definition
% Column 1: Start sample for trials
% Column 2: End sample for trials
% Column 3: Offset of beginning of each trial from Trigger onset in Samples (i.e. -100).
%__________________________________________________________________________
%% === Example 1:
% cfg = [];
% cfg.datafile = 'c,rfDC';% xxx is directory of the raw MEG scan
% cfg.trialfun = 'trialfun_ObjectNaming_Grid';
% cfg.trialdef.stimulCode = [1]; % Object
% cfg.trialdef.prestimTime = 1;
% cfg.trialdef.poststimTime = 1;
% %cfg.trialdef.plotresults = 'yes';%(default = 'no')
% cfgDefTr = ft_definetrial(cfg);
% cfgDefTr.dataformat = '4d';
% cfgDefTr.headerformat = '4d';
% dataRaw = ft_preprocessing(cfgDefTr);

%__________________________________________________________________________

%%
if ~isfield(cfg.trialdef, 'plotresults'),       cfg.trialdef.plotresults = 'no'; end

datafile = cfg.datafile;
stimulCode = [1]; % VY edited 
% stimulCode = cfg.trialdef.stimulCode;
prestimTime = cfg.trialdef.prestimTime;
poststimTime = cfg.trialdef.poststimTime;
if strcmpi(cfg.trialdef.plotresults , 'yes' ), 
    plot_flag = 1;
    titletxt = {'NounVrbPrt','Cross','Resp'};
else
    plot_flag = 0;
end
%%
if stimulCode == 3 && poststimTime < 3%All
    poststimTime = 3.999;
    warning( [ 'The "poststimTime" is not enough to cover Encode , Mask , and Recall periods.' , 13, ' The "poststimTime" is set to ' , num2str(poststimTime) , ' sec.',13 ,13] )
end
%===== trigger code based 
hdr = ft_read_header(datafile);
Fsample = hdr.Fs;
prestimSamples = floor(prestimTime*Fsample);
poststimSamples = floor(poststimTime*Fsample);

TriggerSyncMask = 2^13 - 1;
Base_trigger_code = [ 2048+4 , 6 ];
% % trigger_code = Base_trigger_code( stimulCode );




%==========================================================================
%%
detailedTrigger = ft_read_data(datafile,'chanindx',1,'header',hdr,'eventformat','4d','dataformat','4d');
%%
trg  = bitand( detailedTrigger , TriggerSyncMask );
trg2 = trg;
% Fix problem of photodiode in some subjects
tmp = tokenize(cfg.datafile,'\');
% SubjectID = tmp{end - 5};
SubjectID = tmp{end - 3}; % vy edited
RunID =  tmp{end - 1};

switch  SubjectID
    case  {'C-101' }
%         trg2( 1:(20*hdr.Fs)) = 0;
         trg2 = trg;
        
    case  {'C-121','C-122', 'C-120', 'C-123'}
        Base_trigger_code = [ 2044 , 6 ];
        trg2 = 2048 - trg;
        
    case  {'C-118', 'C-119'}
        Base_trigger_code = [ 4 , 6 ];
        trg2 = trg;
    otherwise
      error('define trigger properly')   
end
trigger_code = Base_trigger_code( stimulCode );

ixp = find( diff(trg2) > 0 ) + 1;
ixp0 = find(trg2( ixp ) == trigger_code );
ixp = ixp( ixp0 );
ixtrlF = ixp(:);


%%
if isempty(ixtrlF)
    disp( [ 'There is no trial in "' cfg.datafile '"'])
    disp( 'based on parameters of following cfg:')
    disp( cfg.trialdef )
    trl = [];
    return;
else
    trl = zeros( length(ixtrlF) , 3);
    trl(:,1) = ixtrlF - prestimSamples;
    trl(:,2) = ixtrlF + poststimSamples;
    trl(:,3) = - prestimSamples;
end

ix0 = find( trl(:,1) <= 0 );
trl(ix0,:) = [];

ix0 = find( trl(:,2) > length(trg) );
trl(ix0,:) = [];

%%
if plot_flag
    try
    t = (0:(hdr.nSamples-1))/Fsample;
    ix = [];
    for i = 1 : size(trl,1) , ix = [ ix, trl(i,1):trl(i,2) ]; end
    trg2 = zeros(length(trg),1);
    trg2(ix) = max(trg(ix));
    h = figure;
    set( h , 'Name' , cfg.datafile );
    plot( t,trg,t,trg2,'r' )
    legend( 'Trigger channel' , 'Selected trigger' )
    title( ['Stimuli: ' titletxt{stimulCode} ] )
    xlabel( 'Time (s)' )
    catch
        warning(['Error in Plot ' titletxt{stimulCode} ])
    end
end

%%
end
