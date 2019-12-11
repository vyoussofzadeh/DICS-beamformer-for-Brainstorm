function    trl = trialfun_KaleidR( cfg )

% This is the Trial Definition Function for Kaleidoscope visual/memory task (ISI ~= 7 sec)
% Each epoch of stimuli contains: 
%       1- Encode:  image for 1.5 sec; Code 2056
%       2- Retain:    image for 1.5 sec; Code 2064
%       3- Recall:  image for 1.5 sec; Code 2080(match) or 2112(not match)
%       4- All ( Encode but poststimTime is long (5s) to cover Retain and Recall) 
%       5- Response(pp):  Code 8(match) or 16(not match)
%
% The trigger is derived from the Photodiode trigger channle
%__________________________________________________________________________
%% Input:
% cfg: Structure with fields that contain parameters for extracting trials
% cfg.datafile: Char string represents filename of raw input MEG data
% cfg.trialdef.stimulCode: Respresting the desired stimulus type to be extracted from the data with following integer values:
% cfg.trialdef.stimulCode = [1,2,3,4, 5];
% % trigger_name = {'Encode','Retain','Recall','All','Response'};
% % trigger_code: see above

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
% cfg.trialfun = 'trialfun_KaleidR';
% cfg.trialdef.stimulCode = [1];
% cfg.trialdef.prestimTime = 0.1;
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
    titletxt = { 'Encode' , 'Retain' , 'Recall'  , 'All' , 'Response'};
else
    plot_flag = 0;
end
%%
if stimulCode == 4 && poststimTime < 4.5%All
    poststimTime = 4.9988;
    warning( [ 'The "poststimTime" is not enough to cover Encode , Retain , and Recall periods.' , 13, ' The "poststimTime" is set to ' , num2str(poststimTime) , ' sec.',13 ,13] )
end
    
%===== trigger code based 
hdr = ft_read_header(datafile);
Fsample = hdr.Fs;
prestimSamples = floor(prestimTime*Fsample);
poststimSamples = floor(poststimTime*Fsample);

TriggerSyncMask = 2^13 - 1;
Base_trigger_code = [   2056    2064    2080    2056    8;
                        0       0       2112    0       16];

if stimulCode == 5
    RespFlag = 1;
    response_code = Base_trigger_code(:,stimulCode);
    response_code( response_code==0 ) = [];   
    trigger_code = Base_trigger_code(:,3);
else
    RespFlag = 0;
    trigger_code = Base_trigger_code(:,stimulCode);
    trigger_code( trigger_code==0 ) = [];
end

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
    resp  = bitand( resp0 , 255 );
    ixp0 = find( diff(resp) > 0 ) + 1;
    ixp = [];
    for j = 1 : length( ixp0 )
        if any( resp( ixp0(j) ) == response_code )
            ixp = [ ixp ; ixp0( j ) ];
        end
    end
    
    %     if length(ixp) < 0.2*length(ixtrlF0)
    %         resp  = bitand( resp0 , 2^7 );
    %         ixp = find( diff(resp) > 0 ) + 1;
    %     end
    
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
        trg  = bitand( trg , TriggerSyncMask );
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
