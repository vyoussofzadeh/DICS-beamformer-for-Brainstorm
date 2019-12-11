savepath = ('speech');
if exist(savepath, 'file') == 0, mkdir(savepath), end

%%
load(['f_',subj,'.mat']);

cfg                         = [];
cfg.dataset                 = datafile;
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = f_data.cfg.trialdef.eventtype;
cfg.trialdef.eventvalue     = f_data.cfg.trialdef.eventvalue; % the value of the stimulus trigger for fully incongruent (FIC).
cfg.trialdef.prestim        = 1; % in seconds
cfg.trialdef.poststim       = 3; % in seconds
cfg = ft_definetrial(cfg);

cfg.channel = {'MISC001'};
cfg.hpfreq = 70;
cfg.demean = 'yes';
speech_data = ft_preprocessing(cfg);

%%
tt=[];
for i=1:length(speech_data.trial)
    
    tmp = speech_data.trial{1,i} - mean(speech_data.trial{1,i});
    
    tmp = detrend(tmp);
%     tmp = tmp - detrend_sdata;
    
    tmp1 = zeros(1,length(tmp));
    tmp1(1,200:end-200) = tmp(1,200:end-200);
    
%     [mx, idx] = max(tmp);
    [thres_buf,env, bin] = envelop_hilbert_modified(abs(tmp1));
    
    [a,b] = find(thres_buf > 0.8.*max(thres_buf));
%     hold on 
% %     plot(b,a,'*')
%     y = ylim; % current y-axis limits
%     plot([b(1) b(1)],[y(1) y(2)])

    [d,initCross,finalCross,nextCross,midRef] =  dutycycle(bin);
    if isempty(initCross)
        idx = find(bin > 0); ipoints = idx(1);
    else
        max_idx = intersect(find(initCross < b(1)),find(finalCross > b(1)));
        if isempty(max_idx)
            max_idx = find(initCross < b(1));
            max_idx = max_idx(end);            
        end
        ipoints = round(initCross(max_idx));
    end

%     
%     [mx,idx] = max(d);
    
    %     [ipoints, residual] = findchangepts(abs(hilbert(thres_buf)));
    %     [ipoints, residual] = findchangepts((thres_buf));
%     [ipoints, residual] = findchangepts(thres_buf);
    
    tt(i) = speech_data.time{1}(ipoints);
    ipoints_all(i) = ipoints;
        ipoints_good(i) = ipoints > 0.5.*mean(ipoints_all);
    
    
    figure, plot(speech_data.time{1}, abs(hilbert(speech_data.trial{1,i}))),
    %             figure, plot(speech_data.time{1}, envelope(speech_data.trial{1,i})),
    hold on,
    scaled = (ipoints * (abs(min(speech_data.time{1,i})) + abs(max(speech_data.time{1}))))/length(speech_data.trial{1,i});
    scaled = scaled - abs(min(speech_data.time{1,i}));
    vline(scaled,'g',['speech onset:', num2str(scaled),'sec']),
    box off;
    set(gca,'color','none');
%                 pause
    
end

thre  = speech_data.time{1}(round(0.5.*mean(ipoints_all)));
L = length(tt);
idx_good = find(ipoints_all >= 0.5.*mean(ipoints_all));
idx_bad  = find(ipoints_all < 0.5.*mean(ipoints_all));

figure,
plot(tt, '*'),
hold on
tt1 = tt; tt1(idx_good) = nan;
plot(tt1, 'r*'),
hline(thre,'b',['mean:', num2str(thre),'sec']),
box off;
set(gca,'color','none');
title('Speech onset'),
ylabel('Time (sec)');
xlabel('Trials');
set(gca,'Xtick', 1:L,'XtickLabel',1:L);
xlim([1 L]);
set(gca,'FontSize',10,'XTickLabelRotation',90);
grid

%%
print('speech/speech','-dpng');
badSpeech = table(idx_bad');
textfile_rej = 'speech/bad_speech';
badSpeech.Properties.VariableNames{'Var1'} = 'bSpeechs';
writetable(badSpeech,textfile_rej,'Delimiter',' ');

%%
[C,ia,ib] = intersect(cln_data.sampleinfo(:,1),speech_data.sampleinfo(idx_good,1));

cfg = [];
cfg.trials = ia;
cln_data = ft_preprocessing(cfg, cln_data);


%%


print('ica/ica2','-dpng');