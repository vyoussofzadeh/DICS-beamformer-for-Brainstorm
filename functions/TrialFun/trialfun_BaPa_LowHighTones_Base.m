function    trl = trialfun_BaPa_LowHighTones_Base( cfg )

% This is the Trial Definition Function for Ba/Pa or Low/High tones.
% The trigger is derived from the trigger channle
%__________________________________________________________________________
%% Input:
% cfg: Structure with fields that contain parameters for extracting trials
% cfg.datafile: Char string represents filename of raw input MEG data
% cfg.trialdef.stimulCode: Respresting the desired stimulus type to be extracted from the data with following integer values:
% cfg.trialdef.stimulCode = [1,2,3, 4,5,6 ]; 
% % trigger_name = {'Ba','Pa','BaPa', 'Lw','Hi','LwHi'};
% % trigger_code:
%     8         -> Low Tone
%     16        -> High Tone
%     64        -> Ba
%     128       -> Pa
%     32        -> Cue

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
% cfg.datafile = 'xxx\c,rfDC';% xxx is directory of the raw MEG scan
% cfg.trialfun = 'trialfun_BaPaShortTone';
% cfg.trialdef.stimulCode = [1]; % Ba
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
stimulCode = cfg.trialdef.stimulCode;% 1 to 12;
prestimTime = cfg.trialdef.prestimTime;
poststimTime = cfg.trialdef.poststimTime;
if strcmpi(cfg.trialdef.plotresults , 'yes' ), 
    plot_flag = 1;
    titletxt = {'Ba','Pa','BaPa', 'Lw','Hi','LwHi'};
else
    plot_flag = 0;
end
%%
%===== trigger code based 
% cfg.trialdef.stimulCode = [1,2,3, 4,5,6 ]; 
% % trigger_name = {'Ba','Pa','BaPa', 'Lw','Hi','LwHi'};
TriggerSyncMask = 2^13 - 1;
Base_trigger_code = [ 64,128, 8,16 ];
trigger_code = [];
switch stimulCode
    case 1 % 'Ba'
        trigger_code = 64;
    case 2 % 'Pa'
        trigger_code = 128;
    case 3 % 'BaPa'
        trigger_code = [ 64 128];
        
    case 4 % 'Lw'
        trigger_code = 8;
    case 5 % 'Hi'
        trigger_code = 16;
    case 6 % 'LwHi'
        trigger_code = [ 8 16 ];
    otherwise
     error( [ 'stimulCode was equal to ', num2str(stimulCode) , ', and it has to be 1,2,3,4,5, or 6' ] );
end
%==========================================================================
%%
hdr = ft_read_header(datafile);
Fsample = hdr.Fs;
detailedTrigger = ft_read_data(datafile,'chanindx',1,'header',hdr,'eventformat','4d','dataformat','4d');
prestimSamples = floor(prestimTime*Fsample);
poststimSamples = floor(poststimTime*Fsample);
%%
trg  = bitand( detailedTrigger , TriggerSyncMask );
ixp = find( diff(trg) > 0 ) + 1;
trgp = trg( ixp );
%%
ixtrlF = [];
for j= 1 : length( trgp )
    if any( trgp( j ) == trigger_code )
        ixtrlF = [ ixtrlF ; ixp( j ) ];
    end
end

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
    t = (0:(hdr.nSamples-1))/Fsample;
    ix = [];
    for i = 1 : size(trl,1) , ix = [ ix, trl(i,1):trl(i,2) ]; end
    trg2 = zeros(size(trg));
    trg2(ix) = max(trg(ix));
    h = figure;
    set( h , 'Name' , cfg.datafile );
    plot( t,trg,t,trg2,'r' )
    legend( 'Trigger channel' , 'Selected trigger' )
    title( ['Stimuli: ' titletxt{stimulCode} ] )
    xlabel( 'Time (s)' )
end

%%
end
