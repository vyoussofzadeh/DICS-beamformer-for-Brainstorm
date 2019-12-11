function h = plot_conn(npsi,label, tit)

% figure, 
h = imagesc(npsi);
set(gca,'YTick', 1:size(npsi,1),'YTickLabel',label);
set(gca,'XTick', 1:size(npsi,1),'XTickLabel',label);
% set(gca,'FontSize',10,'XTickLabelRotation',90)
colorbar
grid on
title(tit);
end