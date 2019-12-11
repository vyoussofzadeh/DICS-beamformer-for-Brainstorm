function    trl = trialfun_VerbGenAud( cfg )

% This is the Trial Definition Function for CRM Auditory Words(x135)experiment.
% The trigger is derived from the trigger channle
%__________________________________________________________________________
%% Input:
% cfg: Structure with fields that contain parameters for extracting trials
% cfg.datafile: Char string represents filename of raw input MEG data
% cfg.trialdef.stimulCode: Respresting the desired stimulus type to be extracted from the data with following integer values:
% % trigger_name = { 'NounVrbAud' };
% % trigger_code = { 64};

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
% cfg.trialfun = 'trialfun_CRM_Grid';
% cfg.trialdef.prestimTime = 1;
% cfg.trialdef.poststimTime = 1;
% %cfg.trialdef.plotresults = 'no';%(default = 'no')
% cfgDefTr = ft_definetrial(cfg);
% cfgDefTr.dataformat = '4d';
% cfgDefTr.headerformat = '4d';
% dataRaw = ft_preprocessing(cfgDefTr);

%__________________________________________________________________________

%%
if ~isfield(cfg.trialdef, 'plotresults'),       cfg.trialdef.plotresults = 'no'; end

datafile = cfg.datafile;
stimulCode = 64;
prestimTime = cfg.trialdef.prestimTime;
poststimTime = cfg.trialdef.poststimTime;
if strcmpi(cfg.trialdef.plotresults , 'yes' ), 
    plot_flag = 1;
    titletxt = { 'NounVrbAud'}
else
    plot_flag = 0;
end
%%
%===== trigger code based 
%Trigger:       Words
%Trigger Code:  64
TriggerSyncMask = 2^13 - 1;
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
    if any( trgp( j ) == stimulCode )
        ixtrlF = [ ixtrlF ; ixp(j) ];
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
%%

ix = find( trl(:,1) <=0 );
ix = [ ix , find( trl(:,2) > length(trg) ) ];
trl(ix,:) = [];

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
    title( ['Stimuli: All Words' ] )
    xlabel( 'Time (s)' )
end

%%
end
