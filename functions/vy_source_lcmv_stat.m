%% leadfield (indv)
% cfg = [];
% cfg.headmodel = individual_headmodel;
% cfg.reducerank = 2;
% cfg.normalize = 'yes';
% cfg.normalizeparam = 0.5;
% cfg.grid.resolution = 10;
% cfg.grid.unit = 'mm';
% individual_grid2 = ft_prepare_leadfield(cfg,t_data.all);
% %
% s_data = vy_source_stat(t_data, individual_grid2, individual_headmodel);
% stats = vy_source_stat_montcarlo(s_data);

%% set threshold
% projthresh = 0.07;
% stats1 = stats;
% val = stats1.stat;
% val(abs(val) < projthresh*max(abs(val(:)))) = NaN;
% stats1.stat = val;
%%
% param = [];
% param.mask = 'stat';
% param.loc = 'min';
% source_int_lcmv = vy_source_plot(stats1,mri_realigned,param,1);
% %
% vy_mapvisualisation(source_int_lcmv,'stat',0, []);
% vy_savenifti(source_int_lcmv,'stat','test.nii');

%% source- single trial
s_data = vy_source_stat(t_data, individual_grid, individual_headmodel);

%%
stats = vy_source_stat_montcarlo(s_data);
stats.pos     = template_grid.pos;
stats.dim     = template_grid.dim;
stats.inside  = template_grid.inside;

stats1  = stats;
stats1.stat =  stats1.stat + 1;

%% Saving data
mtd = 'lcmv-stat';
outputdir_lcmvstat = fullfile(outputdir,mtd);
if exist(outputdir_lcmvstat, 'file') == 0, mkdir(outputdir_lcmvstat), end
savedata = fullfile(outputdir_lcmvstat,['s_lmcv_stat',subj,'_',run,'.mat']);
save(savedata, 'stats1', '-v7.3');

%%
param = [];
param.mask = 'stat';
param.loc = 'min';
source_int_lcmv = vy_source_plot(stats1,template_mri,param,2);
savefig = fullfile(outputdir_lcmvstat,[mtd,'_1_',subj,'_',run]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

clear savepath
savepath{1} = fullfile(outputdir_lcmvstat,[mtd,'_par_2',subj,'_',run]);
savepath{2} = fullfile(outputdir_lcmvstat,[mtd,'_par_3',subj,'_',run]);
vy_mapvisualisation(source_int_lcmv,'stat',0.6, savepath);

savenii = fullfile(outputdir_lcmvstat,['s_lcmv_stat_',subj,'_',run,'.nii']);
vy_savenifti(source_int_lcmv,'stat',savenii);

%% parcellation
% [~, data_intpar, coor] = vy_parcellate(stats, atlas,'stat');
% 
% savepath = fullfile(outputdir_lcmvstat,[mtd,'_par_1',subj,'_',run,'.mat']);
% save(savepath, 'data_intpar', '-v7.3');
% 
% % clear savepath
% % savepath{1} = fullfile(outputdir_lcmvstat,[mtd,'_par_2',subj,'_',run]);
% % savepath{2} = fullfile(outputdir_lcmvstat,[mtd,'_par_3',subj,'_',run]);
% % vy_mapvisualisation(data_intpar,'pow',0.6, savepath)
% vy_mapvisualisation(data_intpar,'stat',0.6, [])
% 
% 
% % save nii
% savenii = fullfile(outputdir_lcmvstat,['s_lcmv_par',subj,'_',run,'.nii']);
% vy_savenifti(data_intpar,'pow',savenii);
% 
% % ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.5, coor, 'stat');
% % disp(ROI_sel)
% savepath = fullfile(outputdir_lcmvstat,['s_dics_ROIs_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

%% 

