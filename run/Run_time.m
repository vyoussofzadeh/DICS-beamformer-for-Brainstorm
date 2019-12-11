%% Epoching
%     switch task
%         case 1 %'DefNam'
%             toi = input('Eneter toi (e.g. [-0.3,0;1.5,2]): ');
%         case 2 %'PicNam'
%             %             toi = [-0.3,0;0.7,1];
%             toi = [-0.2,0;0.5,0.8];
%     end
%     toi = [-0.2,0;1,1.5];
% toi = [-0.3,0;0,1.5];
ep_data = vy_epoch(datain, toi);

%- Appending data
cfg = [];
ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);

%% Data-covarinace estimation
t_data = vy_timelock(ep_data);

%%
% savepath = fullfile(outd.sub,'Timelock',['t_',subj,'.mat']);
% tlk = t_data.all;
% save(savepath, 'tlk', '-v7.3');
