
% freq analysis - prepration for DICS source analysis
f_data.bsl = vy_fft(ep_data.bsl, [2,40], 0,[],0); f_data.bsl.elec = sens;
f_data.pst = vy_fft(ep_data.pst, [2,40], 0,[],0); f_data.pst.elec = sens;

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

outputdir_dics = fullfile(outputdir,'dics');
if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end

savepath = fullfile(outputdir_dics,['psd_',subj,'_',run,'.mat']);
save(savepath, 'ff','psd_bsl','psd_pst', '-v7.3');
savepath = fullfile(outputdir_dics,['fftcon_',subj,'_',run]);
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

% [a,b] = min(psd_pst - psd_bsl);
% f = ff(b); f = round(f);

f = input('Freq of intereset (Hz)? ');
% Freq of interest - prepration for DICS source analysis
f_data.app = vy_fft(ep_data.app, [f,f], 0,[],0); f_data.app.elec = sens;
f_data.bsl = vy_fft(ep_data.bsl, [f,f], 0,[],0); f_data.bsl.elec = sens;
f_data.pst = vy_fft(ep_data.pst, [f,f], 0,[],0); f_data.pst.elec = sens;
%%
[s_data_dics, ~] = vy_source_freq(f_data, individual_grid, individual_headmodel, 'dics');

%%
cfg = [];
cfg.parameter = 'pow';
cfg.operation = '(x1-x2)/(x1+x2)';
source_diff_dics = ft_math(cfg,s_data_dics.pst,s_data_dics.bsl);
source_diff_dics.pos     = template_grid.pos;
source_diff_dics.dim     = template_grid.dim;
source_diff_dics.inside  = template_grid.inside;

%%
% outputdir_dics = fullfile(outputdir,'dics');
% if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end
% savedata = fullfile(outputdir_dics,['s_dics_',subj,'_',run,'.mat']);
% save(outputdir_dics, 'source_diff_dics', '-v7.3');

mtd = 'source_dics';
param = [];
param.mask = 'pow';
param.loc = 'min';
source_int_dics = vy_source_plot(source_diff_dics,template_mri,param,2);
% savefig = fullfile(outputdir_dics,[mtd,'_1_',subj,'_',run]);
% hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

% clear savepath
% savepath{1} = fullfile(outputdir_dics,[mtd,'_2_',subj,'_',run]);
% savepath{2} = fullfile(outputdir_dics,[mtd,'_3_',subj,'_',run]);
vy_mapvisualisation(source_int_dics,'pow',0.6, []);

% save nii
% savenii = fullfile(outputdir_dics,['s_dics_',subj,'_',run,'.nii']);
% vy_savenifti(source_int_dics,'pow',savenii);

%% parcellation
[~, data_intpar, coor] = vy_parcellate(source_diff_dics, atlas,'pow');

% savepath = fullfile(outputdir_dics,[mtd,'_par_1',subj,'_',run,'.mat']);
% save(savepath, 'data_intpar', '-v7.3');
% 
% clear savepath
% savepath{1} = fullfile(outputdir_dics,[mtd,'_par_2',subj,'_',run]);
% savepath{2} = fullfile(outputdir_dics,[mtd,'_par_3',subj,'_',run]);
vy_mapvisualisation(data_intpar,'pow',0.7, [])

% save nii
savenii = fullfile(outputdir_dics,['s_dics_par',subj,'_',run,'.nii']);
vy_savenifti(data_intpar,'pow',savenii);

% ROI summary
[ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, 'pow');
disp(ROI_sel)
savepath = fullfile(outputdir_dics,[mtd,'_par_roi',subj,'_',run]);
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);





