function vy_source_lcmv_beta(cfg_main, cln_data)

cfg = [];
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfreq = 12;
cfg.lpfreq = 28;
fcln_data = ft_preprocessing(cfg, cln_data);
hp = cfg.hpfreq;
lp = cfg.lpfreq;

%%
fep_data = vy_epoch(fcln_data, cfg_main.toi);
cfg = [];
fep_data.app = ft_appenddata(cfg,fep_data.bsl,fep_data.pst);
ft_data = vy_timelock(fep_data);

%%
[s_data_lcmv, ~] = vy_source(ft_data, cfg_main.grid, cfg_main.headmodel);

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'log10(x1/x2)';
source_diff_lcmv = ft_math(cfg,s_data_lcmv.pst,s_data_lcmv.bsl);
% load temp_grid
source_diff_lcmv.pos     = cfg_main.template_grid.pos;
source_diff_lcmv.dim     = cfg_main.template_grid.dim;
source_diff_lcmv.inside  = cfg_main.template_grid.inside;
source_diff_lcmv.pow(source_diff_lcmv.pow>0)=0;


%% plotting
outputdir_lcmv = cfg_main.outputdir;
if exist(outputdir_lcmv, 'file') == 0, mkdir(outputdir_lcmv), end

savefig = fullfile(outputdir_lcmv,[num2str(hp),'_',num2str(lp),'Hz','_1_',cfg_main.subj]);

cfg = [];
cfg.mask = 'pow';
cfg.loc = 'min';
cfg.template = cfg_main.template_mri;
cfg.savefile = savefig;
cfg.volnorm     = 2; % yes: 1
source_int_lcmv = vy_source_plot(cfg, source_diff_lcmv);


% clear savepath
% savepath{1} = fullfile(outputdir_lcmv,[num2str(hp),'_',num2str(lp),'Hz','_2_',cfg_main.subj]);
% savepath{2} = fullfile(outputdir_lcmv,[num2str(hp),'_',num2str(lp),'Hz','_3_',cfg_main.subj]);
% vy_mapvisualisation(source_int_dics,cfg.mask,0.6, savepath);
vy_mapvisualisation(source_int_lcmv,cfg.mask,0.6, []);

%% stats
% stat = vy_source_stat_montcarlo(s_data_lcmv);
% stat.pos     = cfg_main.template_grid.pos;
% stat.dim     = cfg_main.template_grid.dim;
% stat.inside  = cfg_main.template_grid.inside;
% 
% tmp = stat.stat;
% tmp2 = zeros(size(stat.pos,1),1);
% tmp2(stat.inside) = tmp;
% 
% stats1  = stat;
% stats1.stat =  tmp2;
% stats1.mask = stat.inside;
% stats1.stat(stats1.stat>0)=0;


end

% outputdir_lcmv = fullfile(outputdir,'lcmv');
% if exist(outputdir_lcmv, 'file') == 0, mkdir(outputdir_lcmv), end
% savedata = fullfile(outputdir_lcmv,['s_lmcv_',subj,'.mat']);
% % save(outputdir_lcmv, 'source_diff_lcmv', '-v7.3');
%
% mtd = 'source_lcmv';
%
% % savefig = fullfile(outputdir_dics,[num2str(f),'Hz','_1_',cfg_main.subj]);
% savefig = fullfile(outputdir   dics,[,subj]);
%
% cfg = [];
% cfg.mask = 'stat';
% cfg.loc = 'min';
% cfg.template = cfg_main.template_mri;
% cfg.savefile = savefig;
% cfg.volnorm     = 2; % yes: 1
% source_int_dics = vy_source_plot(cfg, stats2);
%
%
% param = [];
% param.mask = 'pow';
% param.loc = 'min';
% source_int_lcmv = vy_source_plot(source_diff_lcmv,template_mri,param,2);
% % source_int_lcmv = vy_source_plot(source_diff_lcmv,mri_realigned,param,1);
% % savefig = fullfile(outputdir_lcmv,[mtd,'_1_',subj,'_',run]);
% % hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
%
% % clear savepath
% % savepath{1} = fullfile(outputdir_lcmv,[mtd,'_2_',subj,'_',run]);
% % savepath{2} = fullfile(outputdir_lcmv,[mltd,'_3_',subj,'_',run]);
% % vy_mapvisualisation(source_int_lcmv,'pow',0.6, savepath);
% vy_mapvisualisation(source_int_lcmv,'pow',0.6, []);

% save nii
% savenii = fullfile(outputdir_lcmv,['s_lcmv_',subj,'_',run,'.nii']);
% vy_savenifti(source_int_lcmv,'pow',savenii);
% vy_savenifti(source_int_lcmv,'pow','test.nii');

%% parcellation
% [~, data_intpar, coor] = vy_parcellate(source_diff_lcmv, atlas,'pow');
%
% savepath = fullfile(outputdir_lcmv,[mtd,'_par_1',subj,'_',run,'.mat']);
% save(savepath, 'data_intpar', '-v7.3');
%
% clear savepath
% savepath{1} = fullfile(outputdir_lcmv,[mtd,'_par_2',subj,'_',run]);
% savepath{2} = fullfile(outputdir_lcmv,[mtd,'_par_3',subj,'_',run]);
% vy_mapvisualisation(data_intpar,'pow',0.6, savepath)
%
% % save nii
% savenii = fullfile(outputdir_lcmv,['s_lcmv_par',subj,'_',run,'.nii']);
% vy_savenifti(data_intpar,'pow',savenii);
%
% % ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, 'pow');
% disp(ROI_sel)
% savepath = fullfile(outputdir_lcmv,['s_dics_ROIs_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);