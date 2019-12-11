
clear dics
k = 1;
for i=1:length(files_sel)-1
    disp(files_sel{i})
    d1 = dir(fullfile(DestDirectory,files_sel{i})); files = {d1.name}'; files_sel1 = files(3:end,1);
    for j=1:length(files_sel1)
        d2 = fullfile(DestDirectory,files_sel{i},files_sel1{j},'dics');
        files = dir(fullfile(d2,'source_dics_par*.mat'));
        if isempty(files) == 0
            load(fullfile(files.folder,files.name));
            dics(k,:) = data_intpar.pow;
            data_dis{k} = fullfile(files.folder,files.name);
            k=k+1;
        end
    end
end

outputdir = fullfile(DestDirectory,'Group ave',mtd);
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end

% data_intpar.pow = (mean(dics,1)./std(dics))';
data_intpar.pow = mean(dics,1)';
data_intpar.pow(data_intpar.pow > 0) = 0;
savepath = fullfile(outputdir,'dics_group.mat');
save(savepath, 'data_intpar', '-v7.3');

%% Source vis - quick inspection
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'pow';
figure
ft_sourceplot(cfg,data_intpar);

savepath = fullfile(outputdir,'s_dics_par1_t_1');
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

savepath = fullfile(outputdir,'s_dics_par_t_1');
save(savepath, 'data_intpar', '-v7.3');
% load(savepath)

%%
% data_intpar1 = data_intpar;
% data_intpar1.eigenvector_cent(data_intpar1.eigenvector_cent < 0.7.*max(data_intpar1.eigenvector_cent)) = NaN; % negative effects
% data_intpar1.eigenvector_centdimord = 'chan';

%% Source vis - template
% data_intpar1.eigenvector_cent = 100.*data_intpar1.pow;
msk = 'pow';
clear savepath
savepath{1} = fullfile(outputdir,'s_dics_par1_t_2');
savepath{2} = fullfile(outputdir,'s_dics_par2_t_3');
vy_mapvisualisation(data_intpar,msk,0.5, savepath);
% vy_mapvisualisation(data_intpar,msk,0.5, []);

savenii = fullfile(outputdir,'s_dics_par.nii');
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
savepath = fullfile(outputdir,'s_dics_ROIs');
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

textfile_rej = fullfile(outputdir,'s_dics_ROI_sel');
writetable(ROI_sel,textfile_rej,'Delimiter',' ')

%% Individuals
outputinddir = fullfile(outputdir,'indiv_s_dics');
if exist(outputinddir, 'file') == 0
    mkdir(outputinddir);   %create the directory
end
for i=1:size(dics,1)
    
    Index = strfind(data_dis{i}, '\');
    d = data_dis{i};
    clear savepath
    savepath{1} = fullfile(outputinddir,['1_',d(Index(end)+1:end-4)]);
    savepath{2} = fullfile(outputinddir,['2_',d(Index(end)+1:end-4)]);
    data_intpar.pow = dics(i,:)';
    data_intpar.pow(data_intpar.pow > 0) = 0;
%     ft_mapvisualisation(msk,data_intpar, [], [], 0,0,savepath,1);
    vy_mapvisualisation(data_intpar,msk,0, savepath);
    
    savenii = fullfile(outputinddir,[d(Index(end)+1:end-4),'.nii']);
    vy_savenifti(data_intpar,msk,savenii);
    
    close all
end
