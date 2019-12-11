
%%
clear pow names data_dis k
pow = [];
k = 1;
for i=1:length(files_sel)
    load(files_sel{i});
    datafile = files_sel{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    Date = datafile(Index(5)+1:Index(6)-1);
    Subj  = datafile(Index(6)+1:Index(7)-1);
    pow(i,:) = source_diff_dics.pow;
    disp(datafile)
    disp(['subj:',Subj])
    disp(['Date:',Date])
    disp([num2str(i),'/', num2str(length(files_sel))])
    disp('============');
    names{i} = [num2str(i),'-',Subj,'-',Date];
    Sub_all{i} = Subj;
end
%-
save([tsk,'_pow.mat'],'pow', 'names','Sub_all','-v7.3');
load([tsk,'_pow.mat'])

%%
clear b
for i=1:length(names)
    b{i} = num2str(i);
end
disp(table([b',names']))

%% outliers
%         [ri, idx] = vy_outlier(evc);
[ri, idx] = vy_outlier_baseline(pow,1);
lout = round(length(idx)/4); % 1/5 of the outliers are thrown
outliers_ID = (names(idx(1:lout)));
disp('bad data:')
disp(outliers_ID');
disp(idx(1:lout)');
Good_ID = (names(idx(end:-1:end-lout)));
disp('Best data:')
disp(Good_ID');
disp(idx(end:-1:end-lout)');
pow_sel = pow(idx(lout+1:end),:);
% pow_sel = pow;
%%
idx1 = idx(lout+1:end);
clear data_dis_sel names_sel
for i=1:size(pow_sel,1)
    data_dis_sel{i} = files_sel{idx1(i)};
    names_sel{i} = names{idx1(i)};
    sub_sel{i} = Sub_all{idx1(i)};
end

% evc_sel = evc;
% data_dis_sel = data_dis;
% names_sel = names;
% sub_sel = sub;

%%
% evc_sel = evc;
% data_dis_sel = data_dis;
% names_sel = names;
% sub_sel = sub;

%% 
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
msk = 'pow';

%%
if exist(fullfile(outputdir,tsk), 'file') == 0
    mkdir(fullfile(outputdir,tsk));   %create the directory
end

%% Individuals
D = source_diff_dics;
for i=1:size(pow,1)
    D.(msk) = pow(i,:)';
    cfg = [];
    cfg.mask = 'pow';
    cfg.loc = 'min';
    cfg.template = template_mri;
    cfg.savefile = [];
    cfg.volnorm     = 2; % yes: 1
    vy_source_plot(cfg, D);
    set(gcf,'Name',names{i}) %select the name you want
%     print(['./', tsk,'/',names{i}],'-depsc');
    print(['./', tsk,'/',names{i}],'-dpng');
    disp([num2str(i),'/', num2str(size(pow,1))])
    disp(names{i})
    pause(2)
    close all
end

%% G-average
evc_n = [];
D = source_diff_dics;
evcn = squeeze(mean(pow_sel,1)); % zero-mean, divide by std, average
% evcn = vy_normalize(pow);

D.(msk) = evcn';

cfg = [];
cfg.mask = 'pow';
cfg.loc = 'min';
cfg.template = template_mri;
cfg.savefile = [];
cfg.volnorm     = 2; % yes: 1
dics_pow = vy_source_plot(cfg, D);

%%
savepath = fullfile(outputdir,'groupave.mat');
cd(outputdir)
save(savepath, 'D','pow', '-v7.3');

%% Source vis - quick inspection
% cfg              = [];
% cfg.method       = 'ortho';
% cfg.funparameter = msk;
% % figure
% ft_sourceplot(cfg,D);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

% cfg = [];
% cfg.subj = 'group';
% cfg.mask = 'pow';
% cfg.thre = 0.6;
% cfg.savepath = savepath;
% vy_mapvisualisation(cfg, D);
% vy_mapvisualisation(D,msk,0.7, []);

savefig = fullfile(outputdir,[tsk,'_dics_group_1']);
cfg = [];
cfg.mask = 'pow';
% cfg.loc = 'min';
cfg.template = template_mri;
cfg.savefile = savefig;
cfg.volnorm     = 2; % yes: 1
D1 = vy_source_plot(cfg, D);
set(gcf,'name','group','numbertitle','off')

%%
clear savepath
savepath{1} = fullfile(outputdir,[tsk,'dics_group_2']);
savepath{2} = fullfile(outputdir,[tsk,'dics_group_3']);

cfg = [];
cfg.subj = 'group';
cfg.mask = 'pow';
cfg.thre = 0.8;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, D1);

%%
% D.eigenvector_cent(D.eigenvector_cent < 0.7.*max(D.eigenvector_cent)) = NaN; % negative effects
% D.eigenvector_centdimord = 'chan';

%% Source vis - template, quick and dirty visualisation
% D.eigenvector_cent = 100.*D.eigenvector_cent;
% gtm = 'eigenvector_cent';
% clear savepath
% 
% savepath{1} = fullfile(outputdir,'n_par1_t_2_wb');
% savepath{2} = fullfile(outputdir,'n_par2_t_3_wb');
savenii = fullfile(outputdir,[tsk,'_groupave.nii']);
% 
% param = [];
% param.mask = 'eigenvector_cent';
% param.loc = 'max';
% D1 = vy_source_plot(D,template_mri,param,2);
% vy_mapvisualisation(D1,'eigenvector_cent',0.6, savepath);
vy_savenifti(D1,msk,savenii);

%% more detailled visualisation
[~, D_par, coor] = vy_parcellate(D, atlas, msk);
D_par.powdimord = 'chan';

savepath = fullfile(outputdir,[tsk,'_group_par.mat']);
save(savepath, 'D_par', '-v7.3');

%%
clear savepath
savepath{1} = fullfile(outputdir,[tsk,'dics_par_group_2']);
savepath{2} = fullfile(outputdir,[tsk,'dics_par_group_3']);
cfg = [];
cfg.subj = [tsk,'_group'];
cfg.mask = 'pow';
cfg.thre = 0.6;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, D_par);

%%
savenii = fullfile(outputdir,[task,'group_par.nii']);
vy_savenifti(D_par, msk, savenii);

[ROI, ROI_sel] = vy_ROI_report(D_par,.8, coor, msk);
% disp(ROI_sel)
% savepath = fullfile(outputdir,'n_ROIs_wb_par');
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

% textfile_rej = fullfile(outputdir,'ROI_sel_wb_par');
% writetable(ROI_sel,textfile_rej,'Delimiter',' ');

%% Parcellation - Individuals
clear par_meg_indv
outputinddir = fullfile(outputdir,tsk);
if exist(outputinddir, 'file') == 0
    mkdir(outputinddir);   %create the directory
end

clear savenii
for k=1:size(pow,1)
    
    Index = strfind(files_sel{k}, '/');
    d = names{k};
    D.pow = pow(k,:)';
    
    %-high-res
    cfg = [];
    cfg.parameter = msk;
    cfg.interpmethod = 'sphere_avg';
    cfg.coordsys     = 'mni';
    D1 = ft_sourceinterpolate(cfg, D, template_mri);
    savenii1{k} = fullfile(outputinddir,[d(Index(1)+2:end-3),'.nii']);
    vy_savenifti(D1,msk,savenii1{k});
    
    %-
    close all
    vol_meg = ft_read_mri(savenii1{k});
    [~, par_meg, ~] = vy_parcellate(vol_meg, atlas,'anatomy');
    par_indv(k,:) = par_meg.anatomy;
    
end
save(['./',tsk, '/par_meg'],'par_indv','par_meg','names','coor')
load(['./',tsk, '/par_meg'])

% load('par_meg');

%-parcellation check
addpath(allpath.spm_path);
par_meg.anatomy = mean(par_indv,1)';
vy_parcellate_plot(par_meg, coor, 'net');

%%
data = [];
data.value = par_indv;
data.label = par_meg.label;
LI_meg = vy_laterality(data);
print(['Lat_',tsk],'-dpng')

%% LI, all subjects, 
subj_meg = names_sel;
subj_meg = names;

figure,bar(LI_meg);
view([90 -90])
L = length(LI_meg);
set(gca,'Xtick', 1:L,'XtickLabel',subj_meg);
% set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gcf, 'Position', [1500   500   500  500]);
grid on
set(gca,'color','none');
axis on
xlim([0,L+1])
xlabel('Subj');
ylabel('Laterality');
set(gcf, 'Position', [1500   500   500  500]);
title(tsk)

%%
% par_crm = [-0.15,0.25;
%             0.11, 0.21;
%             0.24, 0.4;
%             -0.11, -0.21;
%             0.40, 0.51;
%             -0.31, 0.31;
%             0.12, 0.47;
%             0.31, 0.25;
%             0.28, 0.43;
%             0.13, 0.27;
%             -0.24, -0.15
%             0.21, 0.31];
% 
% save('laterality_CRM','par_crm');
% load laterality_CRM
figure,
% subplot 121
hbar = bar(LI_meg);
view([90 -90])
L = length(LI_meg);
for i=1:length(idx)
    S{i} = num2str(idx(i));
end
% set(gca,'Xtick',idx,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
set(gcf, 'Position', [800   500   1000  500]);
title(tsk,'fontsize',16)
% legend([{'fMRI-tValue'};{'MEG'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);
for i=1:L
    S{i} = num2str(i);
end
set(gca,'Xtick',1:L,'XtickLabel',S);

%%
% eeg = D_par.eigenvector_cent;
idx = [12:2:20,80:2:88];
idx2 = [idx, idx-1];

% par_eeg.anatomy = eeg2;
% sourceint_pow = vy_parcellate_plot(par_eeg, coor, 'net');
% bar_input = eeg(idx2)';

bar_input = mean(par_indv,1)'; bar_input = bar_input(idx2);
errorbar_input = std(par_indv)'; errorbar_input = errorbar_input(idx2);
% errorbar_input=[fmri_sd(idx2),eeg_sd(idx2)]';

clear label
for i=1:length(idx2)
    label{i}  = par_meg.label{idx2(i)};
end

label = {'Frontal-Inf-Oper-R' ;
    'Frontal-Inf-Tri-R'  ;
    'Frontal-Inf-Orb-R'  ;
    'Rolandic-Oper-R'    ;
    'Supp-Motor-Area-R'  ;
    'Heschl-R'           ;
    'Temporal-Sup-R'     ;
    'Temporal-Pole-Sup-R';
    'Temporal-Mid-R'     ;
    'Temporal-Pole-Mid-R';
    'Frontal-Inf-Oper-L' ;
    'Frontal-Inf-Tri-L'  ;
    'Frontal-Inf-Orb-L'  ;
    'Rolandic-Oper-L'    ;
    'Supp-Motor-Area-L'  ;
    'Heschl-L'           ;
    'Temporal-Sup-L'     ;
    'Temporal-Pole-Sup-L';
    'Temporal-Mid-L'     ;
    'Temporal-Pole-Mid-L'};

L = length(bar_input);
figure, bar(abs(bar_input), 'BarWidth', 0.9);
% figure
set(gca,'Xtick', 1:L,'XtickLabel',label);
% errorbar_groups(bar_input',errorbar_input', 'bar_names',label,'bar_width',0.75,'errorbar_width',0.5);
set(gca,'FontSize',10,'XTickLabelRotation',90);
% grid
box off
set(gca,'color','none');
set(gcf, 'Position', [500   500   1000   400]);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gca,'FontName','Arial');
% legend(tsk);
savepath = fullfile(outputdir,'n_ROIs');
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

roi = bar_input;
save([tsk,'_ROI.mat'], 'roi', '-v7.3');

%%
% figure,bar(LI_meg);
% view([90 -90])
% L = length(LI_meg);
% set(gca,'Xtick', 1:L,'XtickLabel',names);
% % set(gca,'FontSize',10,'XTickLabelRotation',90)
% set(gcf, 'Position', [1500   500   500  500]);
% grid on
% set(gca,'color','none');
% axis on
% xlim([0,L+1])
% xlabel('Subj');
% ylabel('Laterality');
% set(gcf, 'Position', [1500   500   500  500]);
% title(tsk)

%%
% Colors = [217,217,255;89,89,89]/256;
% 
% DataArray = mean(par_indv,1)';
% DataArray = DataArray(idx2);
% 
% figure
% bar_handle = bar(DataArray,'grouped','BarWidth', 1);
% set(bar_handle(1),'FaceColor',Colors(2,:))
% % set(bar_handle(2),'FaceColor',Colors(2,:))
% L = length(DataArray);
% set(gca,'Xtick', 1:L,'XtickLabel',label); set(gca,'FontSize',10,'XTickLabelRotation',45)
% % set(gca,'Xtick', 1:L,'XtickLabel',([1:10,1:10]));
% 
% grid off
% box off
% set(gca,'color','none');
% % axis on
% xlim([0,L+1])
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% xlabel('ROI','FontSize', 20);
% ylabel('Power ERD ','FontSize', 20);
% set(gcf, 'Position', [1000   500   1000  400]);
% % title('VGA')
% legend({tsk})


%% Correlation between meg-net and fMRI
% switch task
%     case 1
%         load('par_fmri_crm.mat');load('par_meg_crm.mat');
%     case 2
%         load('par_fmri_vga.mat');load('par_meg_vga.mat');
% end
% fmri = par_fmri.anatomy./max(par_fmri.anatomy(:));
% meg = par_meg.anatomy./max(par_meg.anatomy(:));
% par_meg.anatomy = meg;
% save('par_meg','par_meg_indv','par_meg','names_sel','coor');
%

%% Laterality index
% close all
% load laterality_CRM
% figure(1),
% subplot 121
% hbar = bar(par_crm);
% view([90 -90])
% L = length(par_crm);
% for i=1:L
%     S{i} = num2str(i);
% end
% % set(gca,'Xtick',idx,'XtickLabel',S);
% xlim([0,L+1]);
% ylim([-1,1]);
% set(gca,'color','none');
% box off
% xlabel('Subj');
% ylabel('Laterality');
% % set(gcf, 'Position', [800   500   1000  500]);
% title('Word-recognition','fontsize',16)
% % legend([{'fMRI-tValue'};{'MEG'}]);
% % set(gca,'FontName','Arial');
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% ylabel('Laterality','FontSize', 16);
% xlabel('Subj','FontSize', 16);
% 
% DataArray = par_crm;
% Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
% figure(2)
% subplot 121
% [xPositions, yPositions, ~, ~] = UnivarScatter(par_crm,'Label',{'fMRI','meg'},'MarkerFaceColor',Colors);
% ylabel('Laterality','FontSize', 16);
% xlabel('Modality','FontSize', 16);
% set(gca,'color','none');
% title('Word-recognition','fontsize',16)
% disp(['fmri: ',num2str(mean(par_crm(:,1))), '+-', num2str(std(par_crm(:,1)))]);
% disp(['meg: ',num2str(mean(par_crm(:,2))), '+-', num2str(std(par_crm(:,2)))]);
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% hold on
% f = [xPositions, yPositions];
% for j=1:length(f)
%     line([f(j,1),f(j,2)],[f(j,3),f(j,4)]);
% end
% 
% 
% for k = 1:numel(hbar)
%     set(hbar(k),'FaceColor',Colors(k,:))
% end
% 
% %
% load laterality_VGA
% figure(1)
% subplot 122
% hbar = bar(par_vga);
% view([90 -90])
% L = length(par_vga);
% for i=1:length(L)
%     S{i} = num2str(i);
% end
% % set(gca,'Xtick',idx,'XtickLabel',S);
% xlim([0,L+1]);
% ylim([-1,1]);
% set(gca,'color','none');
% box off
% xlabel('Subj');
% ylabel('Laterality');
% % set(gcf, 'Position', [800   500   1000  500]);
% title('Verb-Generation','fontsize',16)
% % legend([{'fmri-tValue'};{'meg-hubs'}]);
% % set(gca,'FontName','Arial');
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% disp(['fmri: ',num2str(mean(par_vga(:,1))), '+-', num2str(std(par_vga(:,1)))]);
% disp(['meg: ',num2str(mean(par_vga(:,2))), '+-', num2str(std(par_vga(:,2)))]);
% % colormap(Colors);
% ylabel('Laterality','FontSize', 16);
% xlabel('Subj','FontSize', 16);
% 
% 
% for k = 1:numel(hbar)
%     set(hbar(k),'FaceColor',Colors(k,:))
% end
% 
% figure(2)
% subplot 122
% [xPositions, yPositions, Label, RangeCut] = UnivarScatter(par_vga,'Label',{'fmri','meg'},'MarkerFaceColor',Colors);
% ylabel('Laterality','FontSize', 16);
% set(gca,'color','none');
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% title('Verb-Generation','fontsize',16)
% xlabel('Modality','FontSize', 16);
% f = [xPositions, yPositions];
% for j=1:length(f)
%     line([f(j,1),f(j,2)],[f(j,3),f(j,4)]);
% end
