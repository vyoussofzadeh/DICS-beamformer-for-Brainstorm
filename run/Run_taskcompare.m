% clear, 
close all,

%%
PN = load('./PN/par_meg.mat');
PN2 = load('PN_pow');
PN3 = load('PN_ROI');

DFN = load('./DFN/par_meg.mat');
DFN2 = load('DFN_pow');
DFN3 = load('DFN_ROI');

%%
[C,ia,ib] = intersect(DFN2.Sub_all,PN2.Sub_all, 'stable');
disp(C')
a = [ia,ib]

par_DFN_sel  = DFN.par_indv(ia,:);
par_PN_sel = PN.par_indv(ib,:);

%% 20 Fronto-temporal regions
idx = [12:2:20,80:2:88];
idx2 = [idx-1, idx];

clear label
for i=1:length(idx2)
    label{i}  = PN.par_meg.label{idx2(i)};
end

label = { 
    'Frontal-Inf-Oper-L' ;
    'Frontal-Inf-Tri-L'  ;
    'Frontal-Inf-Orb-L'  ;
    'Rolandic-Oper-L'    ;
    'Supp-Motor-Area-L'  ;
    'Heschl-L'           ;
    'Temporal-Sup-L'     ;
    'Temporal-Pole-Sup-L';
    'Temporal-Mid-L'     ;
    'Temporal-Pole-Mid-L';
    'Frontal-Inf-Oper-R' ;
    'Frontal-Inf-Tri-R'  ;
    'Frontal-Inf-Orb-R'  ;
    'Rolandic-Oper-R'    ;
    'Supp-Motor-Area-R'  ;
    'Heschl-R'           ;
    'Temporal-Sup-R'     ;
    'Temporal-Pole-Sup-R';
    'Temporal-Mid-R'     ;
    'Temporal-Pole-Mid-R';};

%%
%-PN
data = [];
data.value = abs(par_PN_sel);
data.label = PN.par_meg.label;
LI_PN = vy_laterality(data);

%-DFN
data = [];
data.value = abs(par_DFN_sel);
data.label = DFN.par_meg.label;
LI_DFN = vy_laterality(data);

%% -
clear DataArray
DataArray = [LI_PN,LI_DFN];

Colors = [0.9 0 0;0 0.9 0];
figure(1)
UnivarScatter(DataArray,'Label',{'PN','DFN'},'MarkerFaceColor',Colors);
ylabel('Laterality','FontSize', 20);
set(gca,'color','none');

% figure,
% plotSpread(DataArray,'distributionMarkers',{'o'},'xNames', {'PN','DFN'},'categoryColors',{'r'});
% ylabel('Laterality','FontSize', 20);
% set(gcf, 'Position', [500   500   500  500]);
% set(gca,'color','none');

figure(2),bar(LI_PN);
view([90 -90])
L = length(LI_PN);
set(gca,'Xtick', 1:L,'XtickLabel',C);
% set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gcf, 'Position', [1500   500   500  500]);
grid on
set(gca,'color','none');
axis on
xlim([0,L+1])
xlabel('Subj');
ylabel('Laterality');
set(gcf, 'Position', [1500   500   500  500]);
title('PN')

figure(3),bar(LI_DFN);
view([90 -90])
L = length(LI_DFN);
set(gca,'Xtick', 1:L,'XtickLabel',C);
% set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gcf, 'Position', [1500   500   500  500]);
grid on
set(gca,'color','none');
axis on
xlim([0,L+1])
xlabel('Subj');
ylabel('Laterality');
set(gcf, 'Position', [1500   500   500  500]);
title('DFN')

%% Laterality index 
%-PN vs DFN
% close all; clc
figure(4),
% subplot 121
hbar = bar(DataArray);
view([90 -90])
L = length(DataArray);
for i=1:L
    S{i} = num2str(i);
end
set(gca,'Xtick',1:L,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
% title('Word-recognition','fontsize',16)
legend([{'PN'};{'DFN'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);
set(gcf, 'Position', [500   500  400  500]);

Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
Colors = [0  0  1;0    1    0.9];

figure(5)
% subplot 122
[xPositions, yPositions, Label, RangeCut] = UnivarScatter(DataArray,'Label',{'PN','DFN'},'MarkerFaceColor',Colors);

ylabel('Laterality','FontSize', 16);
xlabel('Modality','FontSize', 16);
set(gca,'color','none');
% title('Word-recognition','fontsize',16)
disp(['PN: ',num2str(mean(DataArray(:,1))), '+-', num2str(std(DataArray(:,1)))]);
disp(['DFN: ',num2str(mean(DataArray(:,2))), '+-', num2str(std(DataArray(:,2)))]);
set(gca,'FontName','HelveticaNeueLT Std Lt');
hold on
f = [xPositions, yPositions];
for j=1:length(f)
   line([f(j,1),f(j,2)],[f(j,3),f(j,4)]); 
end

for k = 1:numel(hbar)
set(hbar(k),'FaceColor',Colors(k,:))
end
set(gcf, 'Position', [500 500 400 500]);

%% outliers
df = diff(DataArray')'; baddata = find(abs(df) > 0.7);
DataArray(baddata,:) = [];

figure(6),
% subplot 121
hbar = bar(DataArray);
view([90 -90])
L = length(DataArray);
for i=1:L
    S{i} = num2str(i);
end
set(gca,'Xtick',1:L,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
% title('Word-recognition','fontsize',16)
legend([{'PN'};{'DFN'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);
set(gcf, 'Position', [500   500  400  500]);

Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
Colors = [0  0  1;0    1    0.9];

figure(7)
% subplot 122
[xPositions, yPositions, Label, RangeCut] = UnivarScatter(DataArray,'Label',{'PN','DFN'},'MarkerFaceColor',Colors);

ylabel('Laterality','FontSize', 16);
xlabel('Modality','FontSize', 16);
set(gca,'color','none');
% title('Word-recognition','fontsize',16)
disp(['PN: ',num2str(mean(DataArray(:,1))), '+-', num2str(std(DataArray(:,1)))]);
disp(['DFN: ',num2str(mean(DataArray(:,2))), '+-', num2str(std(DataArray(:,2)))]);
set(gca,'FontName','HelveticaNeueLT Std Lt');
hold on
f = [xPositions, yPositions];
for j=1:length(f)
   line([f(j,1),f(j,2)],[f(j,3),f(j,4)]); 
end

for k = 1:numel(hbar)
set(hbar(k),'FaceColor',Colors(k,:))
end
set(gcf, 'Position', [500 500 400 500]);

%% Visualise parcellation

PN.par_indv(baddata,:) = [];
addpath(allpath.spm_path);
par_meg.anatomy = mean(PN.par_indv,1)';
vy_parcellate_plot(par_meg, coor, 'PN_par');
PN.roi = par_meg.anatomy(idx2);

DFN.par_indv(baddata,:) = [];
par_meg.anatomy = mean(DFN.par_indv,1)';
vy_parcellate_plot(par_meg, coor, 'DFN_par');
DFN.roi = par_meg.anatomy(idx2);


%% Group
clear DataArray
% DataArray = [par_eeg_indv_CRM_sel(1,idx2)',par_eeg_indv_VGA_sel(1,idx2)'];
DataArray = ([DFN.roi,PN.roi]);
DataArray1 = abs(DataArray);

Colors = [217,217,255;89,89,89]/256;

figure(6)
bar_handle = bar(DataArray1,'grouped','BarWidth', 1);
set(bar_handle(1),'FaceColor',Colors(1,:))
set(bar_handle(2),'FaceColor',Colors(2,:))
L = length(DataArray1);
set(gca,'Xtick', 1:L,'XtickLabel',label); set(gca,'FontSize',10,'XTickLabelRotation',45)
% set(gca,'Xtick', 1:L,'XtickLabel',([1:10,1:10]));

r = corr2(DataArray1(:,1),DataArray1(:,2))

grid off
box off
set(gca,'color','none');
% axis on
xlim([0,L+1])
set(gca,'FontName','HelveticaNeueLT Std Lt');
xlabel('ROI','FontSize', 20);
ylabel('DICS-BF abs(power)','FontSize', 20);
set(gcf, 'Position', [1000   500   1000  400]);
title('Group')
legend({'PN','DFN'})

coor1 = coor(idx2,:);

% [idx, l] = sort(DataArray1(:,1),'descend');
% label(l(1:5))
% 13*idx(1:5);
% coor1(l(1:5),:);
% 
% [idx, l] = sort(DataArray1(:,2),'descend');
% label(l(1:5))
% 13.*idx(1:5);
% coor1(l(1:5),:);

%% Parcellation, individual
subj = 1;
DFN_sub = DFN.par_indv(subj,idx2)';
PN_sub = PN.par_indv(subj,idx2)';

clear DataArray
DataArray = ([DFN_sub,PN_sub]);
DataArray1 = abs(DataArray);

Colors = [217,217,255;89,89,89]/256;

figure
bar_handle = bar(DataArray1,'grouped','BarWidth', 1);
set(bar_handle(1),'FaceColor',Colors(1,:))
set(bar_handle(2),'FaceColor',Colors(2,:))
L = length(DataArray1);
set(gca,'Xtick', 1:L,'XtickLabel',label); set(gca,'FontSize',10,'XTickLabelRotation',45)
grid off
box off
set(gca,'color','none');
% axis on
xlim([0,L+1])
set(gca,'FontName','HelveticaNeueLT Std Lt');
xlabel('ROI','FontSize', 20);
ylabel('DICS-BF abs(power)','FontSize', 20);
set(gcf, 'Position', [1000   500   1000  400]);
title('Subject')
legend({'DFN','PN'})

%-DFN
data = [];
data.value = abs(DFN.par_indv(subj,:)); 
data.value = data.value./max(data.value);
data.label = DFN.par_meg.label;
LI_sub_DFN = vy_laterality(data);

%-PN
data = [];
data.value = abs(PN.par_indv(subj,:));
data.value = data.value./max(data.value);
data.label = PN.par_meg.label;
LI_sub_PN = vy_laterality(data);

DataArray = [LI_sub_DFN,LI_sub_PN];
figure,
bar(DataArray);
L = length(DataArray);
for i=1:L
    S{i} = num2str(i);
end
set(gca,'Xtick',1:L,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
% title('Word-recognition','fontsize',16)
legend([{'DFN'};{'PN'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);
set(gcf, 'Position', [500   500  400  500]);

%%

