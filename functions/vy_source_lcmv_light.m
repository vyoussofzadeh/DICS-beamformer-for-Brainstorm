
cfg = [];
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfreq = 12;
cfg.lpfreq = 28;
fcln_data = ft_preprocessing(cfg, cln_data);
%%
fep_data = vy_epoch(fcln_data, toi);
cfg = [];
fep_data.app = ft_appenddata(cfg,fep_data.bsl,fep_data.pst);
ft_data = vy_timelock(fep_data);

%%

[s_data, s_data2] = vy_source(ft_data, individual_grid, individual_headmodel);

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'log10(x1/x2)';
source_diff_lcmv = ft_math(cfg,s_data.pst,s_data.bsl);
% load temp_grid
source_diff_lcmv.pos     = template_grid.pos;
source_diff_lcmv.dim     = template_grid.dim;
source_diff_lcmv.inside  = template_grid.inside;


%% plotting
outputdir_lcmv = cfg_main.outputdir;
if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end


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