
clear evc conn
k = 1;
for i=1:length(files_sel)
    disp(files_sel{i})
    d1 = dir(fullfile(DestDirectory,files_sel{i})); files = {d1.name}'; files_sel1 = files(3:end,1);
    for j=1:length(files_sel1)
        d2 = fullfile(DestDirectory,files_sel{i},files_sel1{j},'network');
        files = dir(fullfile(d2,'network_evc_C-*.mat'));
        if isempty(files) == 0
            load(fullfile(files.folder,files.name));
            evc(k,:) = network_diff_lcmv.eigenvector_cent;
            data_dis{k} = fullfile(files.folder,files.name);
            k=k+1;
        end
    end
end

outputdir = fullfile(DestDirectory,'Group ave',mtd);
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end
D = network_diff_lcmv;
D.eigenvector_cent = (mean(evc,1)./std(evc))';

savepath = fullfile(outputdir,'eigenvec_group_wb.mat');
save(savepath, 'D', '-v7.3');

%% Source vis - quick inspection
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'eigenvector_cent';
figure
ft_sourceplot(cfg,D);
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
savepath = fullfile(outputdir,'n_par1_t_1_wb');

%%
% D.eigenvector_cent(D.eigenvector_cent < 0.7.*max(D.eigenvector_cent)) = NaN; % negative effects
% D.eigenvector_centdimord = 'chan';

%% Source vis - template
% D.eigenvector_cent = 100.*D.eigenvector_cent;

gtm = 'eigenvector_cent';
clear savepath
savepath{1} = fullfile(outputdir,'n_par1_t_2_wb');
savepath{2} = fullfile(outputdir,'n_par2_t_3_wb');
savenii = fullfile(outputdir,'n_wb.nii');
vy_mapvisualisation(D,gtm,0.5, savepath);
vy_savenifti(D,gtm,savenii);

%%
gtm = 'eigenvector_cent';
[~, D_par, coor] = vy_parcellate(D, atlas, gtm);
D_par.eigenvector_centdimord = 'chan';

savepath = fullfile(outputdir,'eigenvec_group_wb_par.mat');
save(savepath, 'D_par', '-v7.3');

clear savepath
savepath{1} = fullfile(outputdir,'n_par1_t_2_wb_par');
savepath{2} = fullfile(outputdir,'n_par2_t_3_wb_par');
vy_mapvisualisation(D_par,gtm,0.7, savepath);

savenii = fullfile(outputdir,'n_wb_par.nii');
vy_savenifti(D_par, gtm, savenii);

[ROI, ROI_sel] = vy_ROI_report(D_par,.8, coor, gtm);
disp(ROI_sel)
savepath = fullfile(outputdir,'n_ROIs_wb_par');
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

textfile_rej = fullfile(outputdir,'ROI_sel_wb_par');
writetable(ROI_sel,textfile_rej,'Delimiter',' ');


%% Individuals
% switch technique
%     case 1
%         outputinddir = fullfile(outputdir,'indiv_n');
%     case 2
%         outputinddir = fullfile(outputdir,'indiv_n_wb');
% end
% if exist(outputinddir, 'file') == 0
%     mkdir(outputinddir);   %create the directory
% end
% for i=1:size(evc,1)
%
%     Index = strfind(data_dis{i}, '\');
%     d = data_dis{i};
%     clear savepath
%     savepath{1} = fullfile(outputinddir,['1_',d(Index(end)+1:end-4)]);
%     savepath{2} = fullfile(outputinddir,['2_',d(Index(end)+1:end-4)]);
%     D.eigenvector_cent = evc(i,:)';
%     vy_mapvisualisation(D,gtm,0, savepath);
%
%     savenii = fullfile(outputinddir,[d(Index(end)+1:end-4),'.nii']);
%     vy_savenifti(D,gtm,savenii);
%
%     close all
% end

