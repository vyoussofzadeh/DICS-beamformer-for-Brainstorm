function    trl = trialfun_VCRM( cfg )

% This is the Trial Definition Function for VCRM Visual Words experiment (ISI ~= 3.75 sec)
% Stimuli contains: 
%       1- Initial Block: 15 words (to be memorized); Code 2056
%       2- Six Blocks: 15 Memorized Words and 10 new words; 
%                 All Words Code 2176
% The trigger is derived from the Photodiode trigger channle
%__________________________________________________________________________
%% Input:
% cfg: Structure with fields that contain parameters for extracting trials
% cfg.datafile: Char string represents filename of raw input MEG data
% cfg.trialdef.stimulCode: Respresting the desired stimulus type to be extracted from the data with following integer values:
% cfg.trialdef.stimulCode = 1;
% % trigger_name = { 'VCRM' };
% % trigger_code = { 2176 };

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
% cfg.trialfun = 'trialfun_VCRM';
% cfg.trialdef.stimulCode = [1]; % AllWords
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
stimulCode = cfg.trialdef.stimulCode;
prestimTime = cfg.trialdef.prestimTime;
poststimTime = cfg.trialdef.poststimTime;
if strcmpi(cfg.trialdef.plotresults , 'yes' ), 
    plot_flag = 1;
    titletxt = { 'VCRM' };%{ 'VCRM' , 'MemorizedWords' , 'NewWords'  , 'Response'};
else
    plot_flag = 0;
end
%%
%===== trigger code based 
hdr = ft_read_header(datafile);
Fsample = hdr.Fs;
prestimSamples = floor(prestimTime*Fsample);
poststimSamples = floor(poststimTime*Fsample);

TriggerSyncMask = 2^13 - 1;
if stimulCode == 4 
    RespFlag = 1; 
    stimulCode = 1;
else
    RespFlag = 0;
end


% Base_trigger_code = [ 2060 , 2080 ];
trigger_code = 2176;
%==========================================================================
%%
detailedTrigger = ft_read_data(datafile,'chanindx',1,'header',hdr,'eventformat','4d','dataformat','4d');
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

if RespFlag
    %-------- find time of all words
    ixtrlF0 = ixtrlF;
    ixtrlF = [];
    %==========================================================================
    resp0 = ft_read_data(datafile,'chanindx',2,'header',hdr,'eventformat','4d','dataformat','4d');
    resp  = bitand( resp0 , 2^3 );
    ixp = find( diff(resp) > 0 ) + 1;
    
    if length(ixp) < 0.2*length(ixtrlF0)
        resp  = bitand( resp0 , 2^7 );
        ixp = find( diff(resp) > 0 ) + 1;
    end
    
    % Delete noisiance response
    ixp2 = [];
    for i = 1 : length(ixp)
        close_trg_ix = max(find( ixp(i) > ixtrlF0 ));
        if( ~isempty( close_trg_ix ) && ( ((ixp(i) - ixtrlF0(close_trg_ix) )/hdr.Fs ) < 2.5 ) )
            ixp2 = [ ixp2, ixp(i)];
        end
    end
    ixp = ixp2;
    
    %---------------------------
    
    ixtrlF = ixp(:);
    if plot_flag
        trg = ft_read_data(datafile,'chanindx',1,'header',hdr,'eventformat','4d','dataformat','4d');
        trg  = bitand( trg , 511 );
        trg = [ trg ; resp ];
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
    trg2 = zeros(length(trg),1);
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
