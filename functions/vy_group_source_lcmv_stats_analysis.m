
clear vol;
k = 1;
for i=1:length(files_sel)-1
    disp(files_sel{i})
    d1 = dir(fullfile(DestDirectory,files_sel{i})); files = {d1.name}'; files_sel1 = files(3:end,1);
    for j=1:length(files_sel1)
        d2 = fullfile(DestDirectory,files_sel{i},files_sel1{j},'lcmv-stat');
        files = dir(fullfile(d2,'s_lcmv*.nii'));
        if isempty(files) == 0
            s_vol = ft_read_mri(fullfile(files.folder,files.name));
            vol(k,:,:,:) = s_vol.anatomy;
            k=k+1;
        end
    end
end


%%
switch task
    case 1
        vol_name = '.\outputs\meg_source_stat_CRM.nii';
        projthresh = 0.4;
    case 2
        vol_name = '.\outputs\meg_source_stat_VGA.nii';
        projthresh = 0.25;
    case 3
        vol_name = '.\outputs\meg_source_stat_VGP.nii';
        projthresh = 0.6;
end

s_vol.anatomy = vy_normalize(vol);
s_vol = vy_vol_thresh(s_vol,projthresh); % abs
% s_vol1 = vy_vol_thresh_posneg(s_vol,.001); % keep positive effects only!

cfg = [];
cfg.filetype  = 'nifti';
cfg.datatype   = 'uint8'; %'float';
cfg.parameter = 'anatomy';
cfg.filename  = vol_name;
ft_volumewrite(cfg, s_vol)
%%
spmpath = 'F:\My Matlab\SPM\spm12_4\spm12\spm12';
connpath = 'F:\My Matlab\Connectivity\Conn\conn';
addpath(connpath);
addpath(genpath(spmpath));

conn_mesh_display(vol_name, '');
title(['\fontsize{16} {\color{magenta}fmri','}'],'interpreter','tex')

%%
outputdir = fullfile(DestDirectory,'Group ave',mtd);
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end
% data_intpar.pow = (mean(lcmv,1)./std(lcmv))';
data_intpar.pow = mean(lcmv,1)';
savepath = fullfile(outputdir,'lcmv_group.mat');
save(savepath, 'data_intpar', '-v7.3');

%% Source vis - quick inspection
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'pow';
figure
ft_sourceplot(cfg,data_intpar);

savepath = fullfile(outputdir,'s_lcmv_par1_t_1');
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

savepath = fullfile(outputdir,'s_lcmv_par_t_1');
save(savepath, 'data_intpar', '-v7.3');

%% Source vis - template

msk = 'pow';
clear savepath
savepath{1} = fullfile(outputdir,'s_lcmv_par1_t_2');
savepath{2} = fullfile(outputdir,'s_lcmv_par2_t_3');
vy_mapvisualisation(data_intpar,msk,0.5, savepath);

savenii = fullfile(outputdir,'s_lcmv_par.nii');
vy_savenifti(data_intpar,msk,savenii);

%% MNI coordinates
coor = [];
pp = data_intpar.brainordinate;
for i=1:length(pp.tissuelabel)
    idx = pp.tissue == i;
    M = mean(pp.pos(idx,:),1);
    M = round(M*100)/100; % rounding into 2 decimals
    coor(i,:) = round(M);
end
%%
[ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, msk);
disp(ROI_sel)
savepath = fullfile(outputdir,'s_lcmv_ROIs');
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

textfile_rej = fullfile(outputdir,'s_lcmv_ROI_sel');
writetable(ROI_sel,textfile_rej,'Delimiter',' ')

%% Individuals
outputinddir = fullfile(outputdir,'indiv_s_lcmv');
if exist(outputinddir, 'file') == 0, mkdir(outputinddir), end
for i=1:size(lcmv,1)
    
    Index = strfind(data_dis{i}, '\');
    d = data_dis{i};
    clear savepath
    savepath{1} = fullfile(outputinddir,['1_',d(Index(end)+1:end)]);
    savepath{2} = fullfile(outputinddir,['2_',d(Index(end)+1:end)]);
    data_intpar.pow = lcmv(i,:)';
    vy_mapvisualisation(data_intpar,msk,0, savepath);
    
    savenii = fullfile(outputinddir,[d(Index(end)+1:end-4),'.nii']);
    vy_savenifti(data_intpar,msk,savenii);
    
    close all
end
