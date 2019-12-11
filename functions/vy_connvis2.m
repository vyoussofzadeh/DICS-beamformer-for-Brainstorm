function [tedge,ROI] = vy_connvis2(conn_ratio,headmodel,conn_par)

% source_conn_bsl.plvspctrm = source_conn_bsl.plvspctrm(source_conn_bsl.inside,source_conn_bsl.inside);
% source_conn_pst.plvspctrm = source_conn_pst.plvspctrm(source_conn_pst.inside,source_conn_pst.inside);
% source_conn_bsl.plvspctrm(isnan(source_conn_bsl.plvspctrm))=0;
% source_conn_pst.plvspctrm(isnan(source_conn_pst.plvspctrm))=0;
% 
% %%
% cfg           = [];
% cfg.operation = 'x1-x2';
% cfg.parameter = conn_par.idx;
% conn_ratio  = ft_math(cfg, source_conn_pst, source_conn_bsl);
% conn_ratio.pos = source_conn_pst.pos;

aedge =  conn_ratio.plvspctrm;
ROI  = conn_ratio.pos(conn_ratio.inside,:);

aedge(isnan(aedge))=0;
tedge = (aedge.* double(aedge > conn_par.conn_thre.*max(aedge(:))));
% figure, imagesc(tedge), colorbar

figure,
ft_plot_vol(headmodel, 'facecolor', [0,0,0], 'edgecolor', 'none');
alpha 0.1;
view ([-196 56])
hold on
% k = 1;
for i = 1:length(tedge)
    for j = 1:length(tedge)
        if tedge(i,j)> 0
            p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
            p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
            pts = [p1; p2];
            line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
%             conn_idx(k,:) = [i,j]; k=k+1;
%             roi_sel(k,:) = ROI
        end
    end
end
view([90 0])
