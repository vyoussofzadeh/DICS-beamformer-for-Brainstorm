function [ROI, ROI_sel] = vy_ROI_report(data,thre,coor, mask)

m = ft_getopt(data, mask);
label = data.label';
L = length(m);

figure,
b = bar(m);
set(b(1,1),'faceColor',[0.7 0.7 0.7]);
% set(b,'BarWidth',1);
% hold on
% errorbar(m_gtm, var_gtm,'r','linestyle', 'none');
% view([90 -90])
% set(gca,'Xtick', 1:L,'XtickLabel',1:L);
set(gca,'Xtick', 1:L,'XtickLabel',label);
box off
set(gca,'color','none');
xlim([0,L+1])
xlabel('ROI (=116)');
ylabel(mask);
set(gcf, 'Position', [100   100   1500   500]);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gca,'FontName','HelveticaNeueLT Std Lt');
grid
hold on
%% sort
% [val,idx] = sort((m),'descend');

if abs(min(m)) > abs(max(m))
    m1 = -m;
    idx2 = find(m1 > thre.*max(m1));
else
    idx2 = find(m > thre.*max(m));
end
n = length(idx2);
col = ones(n,3);
col(:,[1,3]) = 0;
for i=1:n
    m_sig = zeros(size(m,1),size(m,2));
    m_sig(idx2(i)) = m(idx2(i));
    c = bar(m_sig);
    set(c,'faceColor',col(i,:));
    set(c,'BarWidth',1);
end
title(['ROIs (',num2str(100*thre),'% threshold)']);

%% roi summary
L = length(m);
id = table([1:L]','VariableNames',{'ID'});
Z = table(m,'VariableNames',{mask});
ROI = [id,cell2table(label),Z];
coor_var = table(coor,'VariableNames',{'MNI'});
ROI(:,4) = coor_var;
ROI.Properties.VariableNames{'Var4'} = 'MNI';
%

%% roi selected
[~,idx] = sort((abs(m(idx2))),'descend');
% disp(ROI)
ROI_sel = ROI(idx2(idx),:);
disp(ROI_sel)