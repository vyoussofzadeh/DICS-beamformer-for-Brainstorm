function vy_source_dics_stats(cfg_main, ep_data)

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 40];
cfg.plotflag  = 2;
cfg.tapsmofrq = 1;
cfg.taper    = 'hanning';
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;

% % freq analysis - prepration for DICS source analysis
% f_data.bsl = vy_fft(ep_data.bsl, [2,40], 0,[],0); f_data.bsl.elec = sens;
% f_data.pst = vy_fft(ep_data.pst, [2,40], 0,[],0); f_data.pst.elec = sens;

% PSD - sensor space
psd_bsl = squeeze(mean(mean(abs(f_data.bsl.fourierspctrm),2),1));
psd_pst = squeeze(mean(mean(abs(f_data.pst.fourierspctrm),2),1));
ff = linspace(1, 40, length(psd_pst));

figure,plot(ff,psd_bsl)
hold on
% plot(ff,psd,'g')
plot(ff,psd_pst,'r')
hold on
plot(ff,psd_pst - psd_bsl,'k')
xlabel('Hz'); ylabel('psd'),legend({'bsl','pst','diff'})

% outputdir_dics = fullfile(cfg_main.outputdir,'dics');
outputdir_dics = cfg_main.outputdir;
if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end

% savepath = fullfile(outputdir_dics,['psd_',subj,'_',run,'.mat']);
% save(savepath, 'ff','psd_bsl','psd_pst', '-v7.3');
% savepath = fullfile(outputdir_dics,['fftcon_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

[a,b] = min(psd_pst - psd_bsl);
% f = ff(b); f = round(f);

f = input('FOI? ');
% f = 20;
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [f f];
cfg.plotflag  = 2;
cfg.tapsmofrq = 4;
cfg.taper    = 'dpss';
% cfg.taper    = 'hanning';
[f_data.app,tapsmofrq] = vy_fft(cfg, ep_data.app); f_data.app.elec = cfg_main.sens;
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;

% f = 20;
% % Freq of interest - prepration for DICS source analysis
% f_data.app = vy_fft(ep_data.app, [f,f], 0,[],0); f_data.app.elec = sens;
% f_data.bsl = vy_fft(ep_data.bsl, [f,f], 0,[],0); f_data.bsl.elec = sens;
% f_data.pst = vy_fft(ep_data.pst, [f,f], 0,[],0); f_data.pst.elec = sens;

%%
cfg = [];
cfg.headmodel = cfg_main.headmodel;
% cfg.sourcemodel = cfg_main.sourcemodel;
cfg.grid = cfg_main.grid;
cfg.mtag = cfg_main.mtag;
s_data_dics = vy_source_freq(cfg, f_data);
% [s_data_dics, ~] = vy_source_freq(f_data, cfg_main.grid, cfg_main.headmodel, 'dics_stat');
stat = vy_source_stat_montcarlo(s_data_dics);
stat.pos     = cfg_main.template_grid.pos;
stat.dim     = cfg_main.template_grid.dim;
stat.inside  = cfg_main.template_grid.inside;

tmp = stat.stat;
tmp2 = zeros(size(stat.pos,1),1);
tmp2(stat.inside) = tmp;

stats1  = stat;
stats1.stat =  tmp2;
stats1.mask = stat.inside;

% disp('1: Postive effects')
% disp('2: Negative effects')
% disp('3: Both effects')
% effect = input(':');
% 
stats2 = stats1;
% if effect == 1
%     stats2.stat(stats1.stat<0)=0;
% elseif effect == 2
%     stats2.stat(stats1.stat>0)=0;
% elseif effect == 3
%     stats2.stat = stats1.stat;
% end
stats2.stat(stats2.stat>0)=0;
stats2.stat(isnan(stats2.stat))=0;

%%
% param = [];
% param.mask = 'stat';
% param.loc = 'min';
% source_diff_dics = vy_source_plot(stats1,cfg_main.template_mri,param,2);

%%
% cfg = [];
% cfg.parameter = 'pow';
% cfg.operation = '(x1-x2)/(x1+x2)';
% source_diff_dics = ft_math(cfg,s_data_dics.pst,s_data_dics.bsl);
% source_diff_dics.pos     = cfg_main.template_grid.pos;
% source_diff_dics.dim     = cfg_main.template_grid.dim;
% source_diff_dics.inside  = cfg_main.template_grid.inside;
% source_diff_dics.pow(source_diff_dics.pow>0)=0;

%%
% outputdir_dics = fullfile(outputdir,'dics');
% if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end
% savedata = fullfile(outputdir_dics,['s_dics_',subj,'_',run,'.mat']);
% save(outputdir_dics, 'source_diff_dics', '-v7.3');

% mtd = 'source_dics_stats';
% savefig = fullfile(outputdir_dics,[num2str(f),'Hz','_1_',cfg_main.subj]);
% savefig = fullfile(outputdir_dics,[num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);
savefig = fullfile(outputdir_dics,[num2str(cfg_main.toi(2,1)),'_',num2str(cfg_main.toi(2,2)),'sec_',num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);


cfg = [];
cfg.mask = 'stat';
cfg.loc  = 'min';
cfg.template = cfg_main.template_mri;
cfg.savefile = savefig;
cfg.volnorm  = 2; % yes: 1
source_int_dics = vy_source_plot(cfg, stats2);
pause(1),

% param = [];
% param.mask = 'stat';
% param.loc = 'min';
% source_int_dics = vy_source_plot(stats1,cfg_main.template_mri,param,2);
% savefig = fullfile(outputdir_dics,[num2str(f),'Hz','_1_',cfg_main.subj]);
% % hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
% print(savefig,'-depsc')

clear savepath
savepath{1} = fullfile(outputdir_dics,[num2str(cfg_main.toi(2,1)),'_',num2str(cfg_main.toi(2,2)),'sec_',num2str(f),'Hz','_2_',cfg_main.subj]);
savepath{2} = fullfile(outputdir_dics,[num2str(cfg_main.toi(2,1)),'_',num2str(cfg_main.toi(2,2)),'sec_',num2str(f),'Hz','_3_',cfg_main.subj]);
vy_mapvisualisation(source_int_dics,cfg.mask,0.6, savepath);
% vy_mapvisualisation(source_int_dics,cfg.mask,0.6, []);

% restoredefaultpath
% addpath((cfg_main.allpath.ft_path));
% ft_defaults
% addpath(genpath(cfg_main.allpath.hcp_path));
% addpath(genpath(cfg_main.allpath.cd_org));
% addpath(genpath(cfg_main.allpath.exfig_path));
%
% cfg = [];
% cfg.maskparam = 'pow';
% cfg.save.savepath =  savepath;
% % cfg.saveformat = '-eps';
% cfg.save.saveformat = '-png';
% cfg.save.pixdim     = 12;
% cfg.projthresh      = 0.6;
% vy_surfmap(cfg, source_int_dics);




