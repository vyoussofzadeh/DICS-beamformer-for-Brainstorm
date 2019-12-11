
clear;clc
close all
%{
fl = rdir( 'E:\Data\VNS\MEG\*\*\VerbGenVis\**\c,rfDC');
for i = 13: 13
    fn = fl(i).name;
    cfg = [];
    cfg.datafile = fn;%'E:\Data\VNS\MEG\C-110\VNS003\VerbGenVis\11%18%17@15_38\1\c,rfDC';%'E:\Data\VNS\MEG\C-101\VNS012\VerbGenVis\02%05%18@10_03\2\c,rfDC';% xxx is directory of the raw MEG scan
    cfg.trialfun = 'trialfun_VergGenVis';
    cfg.trialdef.stimulCode = [1]; % Object
    cfg.trialdef.prestimTime = 0;
    cfg.trialdef.poststimTime = .01;
    cfg.trialdef.plotresults = 'yes';%(default = 'no')
    cfgDefTr = ft_definetrial(cfg);
    set(gcf,'name' , fn)
    pause (0.0001)
    % cfgDefTr.dataformat = '4d';
    % cfgDefTr.headerformat = '4d';
    % dataRaw = ft_preprocessing(cfgDefTr);
end

%}

fl = rdir( 'E:\Data\VNS\MEG\*\*\VerbGenAud\**\c,rfDC');
for i = 1: numel(fl)
    fn = fl(i).name;
    cfg = [];
    cfg.datafile = fn;%'E:\Data\VNS\MEG\C-110\VNS003\VerbGenVis\11%18%17@15_38\1\c,rfDC';%'E:\Data\VNS\MEG\C-101\VNS012\VerbGenVis\02%05%18@10_03\2\c,rfDC';% xxx is directory of the raw MEG scan
    cfg.trialfun = 'trialfun_VerbGenAud';
    cfg.trialdef.stimulCode = [1]; % Object
    cfg.trialdef.prestimTime = 0;
    cfg.trialdef.poststimTime = .01;
    cfg.trialdef.plotresults = 'yes';%(default = 'no')
    cfgDefTr = ft_definetrial(cfg);
    set(gcf,'name' , fn)
    pause (0.0001)
    % cfgDefTr.dataformat = '4d';
    % cfgDefTr.headerformat = '4d';
    % dataRaw = ft_preprocessing(cfgDefTr);
    disp(fn)
    disp( size(cfgDefTr.trl,1))
    disp(13)
end