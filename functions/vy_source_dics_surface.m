function [source_diff_dics] = vy_source_dics_surface(cfg_main, ep_data)

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 40];
cfg.plotflag  = 2;
cfg.tapsmofrq       = 1;
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

outputdir_dics = fullfile(cfg_main.outputdir,'dics');
if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end

% savepath = fullfile(outputdir_dics,['psd_',subj,'_',run,'.mat']);
% save(savepath, 'ff','psds_bsl','psd_pst', '-v7.3');
% savepath = fullfile(outputdir_dics,['fftcon_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

[a,b] = min(psd_pst - psd_bsl);
% f = ff(b); f = round(f);

f = input('FOI? '); tapsmofrq = 3;
% f = 20;
% cfg = [];
% cfg.savefile = [];
% cfg.saveflag = 2;
% cfg.foilim = [f f];
% cfg.plotflag  = 2;
% f_data.app = vy_fft(cfg, ep_data.app); f_data.app.elec = cfg_main.sens;
% f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
% f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;


cfg.taper    = 'dpss'; cfg.tapsmofrq  = tapsmofrq;
if f < 4, cfg.tapsmofrq  = 1; cfg.taper    = 'hanning'; end
[f_data.app,~,~,tapsmofrq] = vy_fft(cfg, ep_data.app); f_data.app.elec = cfg_main.sens;
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;


% f = 20;
% % Freq of interest - prepration for DICS source analysis
% f_data.app = vy_fft(ep_data.app, [f,f], 0,[],0); f_data.app.elec = sens;
% f_data.bsl = vy_fft(ep_data.bsl, [f,f], 0,[],0); f_data.bsl.elec = sens;
% f_data.pst = vy_fft(ep_data.pst, [f,f], 0,[],0); f_data.pst.elec = sens;
%%
% [s_data_dics, ~] = vy_source_freq(f_data, cfg_main.grid, cfg_main.headmodel, 'dics');

%%
cfg = [];
cfg.headmodel = cfg_main.headmodel;
cfg.sourcemodel = cfg_main.sourcmodel;
cfg.grid = cfg_main.grid;
cfg.mtag = 'dics_fs';
s_data_dics = vy_source_freq(cfg, f_data);

%%
cfg = [];
cfg.parameter = 'pow';
cfg.operation = '(x1-x2)/(x1+x2)';
% cfg.operation = '(x1-x2)';
source_diff_dics = ft_math(cfg,s_data_dics.pst,s_data_dics.bsl);
% source_diff_dics.pos     = cfg_main.template_grid.pos;
% source_diff_dics.dim     = cfg_main.template_grid.dim;
% source_diff_dics.inside  = cfg_main.template_grid.inside;
source_diff_dics.pow(source_diff_dics.pow>0)=0;

%%
% source_diff_dics.tri = sourcemodelT.tri; %ft_sourceplot need the triangulation information, add to source reconstruction.
% 
% 
% cfg = [];
% cfg.method          = 'surface';
% cfg.funparameter    = 'pow';
% cfg.funcolormap     = 'jet';
% % cfg.latency         = [0.7 1];     % The time-point to plot
% cfg.colorbar        = 'no';
% % cfg.avgovertime     = 'yes';
% ft_sourceplot(cfg, source_diff_dics);
% view([-100,20]);
% 
% %%
% % outputdir_dics = fullfile(outputdir,'dics');
% % if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end
% % savedata = fullfile(outputdir_dics,['s_dics_',subj,'_',run,'.mat']);
% % save(outputdir_dics, 'source_diff_dics', '-v7.3');
% 
% mtd = 'source_dics';
% param = [];
% param.mask = 'pow';
% param.loc = 'min';
% source_int_dics = vy_source_plot(source_diff_dics,cfg_main.template_mri,param,2);
% % savefig = fullfile(outputdir_dics,[mtd,'_1_',subj,'_',run]);
% % hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
% 
% % clear savepath
% % savepath{1} = fullfile(outputdir_dics,[mtd,'_2_',subj,'_',run]);
% % savepath{2} = fullfile(outputdir_dics,[mtd,'_3_',subj,'_',run]);
% % vy_mapvisualisation(source_int_dics,'pow',0.6, savepath);
% vy_mapvisualisation(source_int_dics,'pow',0.6, []);

% save nii
% savenii = fullfile(outputdir_dics,['s_dics_',subj,'_',run,'.nii']);
% vy_savenifti(source_int_dics,'pow',savenii);

%% parcellation
% [~, data_intpar, coor] = vy_parcellate(source_diff_dics, atlas,'pow');
% 
% savepath = fullfile(outputdir_dics,[mtd,'_par_1',subj,'_',run,'.mat']);
% save(savepath, 'data_intpar', '-v7.3');
% 
% clear savepath
% savepath{1} = fullfile(outputdir_dics,[mtd,'_par_2',subj,'_',run]);
% savepath{2} = fullfile(outputdir_dics,[mtd,'_par_3',subj,'_',run]);
% vy_mapvisualisation(data_intpar,'pow',0.6, savepath)
% 
% % save nii
% savenii = fullfile(outputdir_dics,['s_dics_par',subj,'_',run,'.nii']);
% vy_savenifti(data_intpar,'pow',savenii);
% 
% % ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, 'pow');
% disp(ROI_sel)
% savepath = fullfile(outputdir_dics,[mtd,'_par_roi',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);





