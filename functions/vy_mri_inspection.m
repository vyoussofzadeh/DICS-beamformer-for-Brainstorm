% function vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir, saveflag)
function vy_mri_inspection(cfg, ~)


switch cfg.mtd
    case 'vol'
        %%
        figure;
        ft_plot_vol(cfg.headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
        hold on;
        ft_plot_headshape(cfg.headshape);
        ft_plot_mesh(cfg.leadfield.pos(cfg.leadfield.inside, :));
        view ([0 90])
        if cfg.saveflag==1
            savepath = fullfile(cfg.outputmridir,'headshape');
            hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
        end
        
    case 'surf'
        %%
        figure; hold on;
        ft_plot_vol(cfg.headmodel, 'facecolor', 'none'); alpha 0.5;
        ft_plot_mesh(cfg.sourcemodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
        
        %%
        figure; hold on;
        ft_plot_mesh(cfg.sourcemodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
        ft_plot_mesh(cfg.leadfield.pos(cfg.leadfield.inside, :));
        
end
%%
% ft_determine_coordsys(individual_headmodel, 'interactive', 'no')
% hold on; % add the subsequent objects to the same figure
% ft_plot_headshape(headshape);
% view ([-10 40 10]);
% if isempty(saveflag)~=2
%     savepath = fullfile(outputmridir,'headshape2');
%     hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
% end

%% headmodel - inspection
% figure;
% ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% hold on;
% ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
% view ([-10 40 10]);
% if isempty(saveflag)~=2
%     savepath = fullfile(outputmridir,'grid');
%     hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
% end

%%
% ft_determine_coordsys(mri_realigned, 'interactive', 'no')
% hold on;
% ft_plot_headshape(headshape);
% view ([50 80 10])
% ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
% if isempty(saveflag)~=2
%     savepath = fullfile(outputmridir,'gridmri');
%     hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
% end
% disp('---------------')
% disp(['figures were saved at,',outputmridir])

%%
% h = figure;
% subplot('position',[0.01 0.51 0.48 0.48]);hold on;
% ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
% % ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
% ft_plot_headshape(headshape,'vertexsize',5); view(180,-90);
% plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
% subplot('position',[0.51 0.51 0.48 0.48]);hold on;
% % ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
% ft_plot_headshape(headshape,'vertexsize',5); view(0,90);
% plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
% subplot('position',[0.01 0.01 0.48 0.48]);hold on;
% ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
% ft_plot_headshape(headshape,'vertexsize',5); view(90,0);
% plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
% subplot('position',[0.51 0.01 0.48 0.48]);hold on;
% % ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
% ft_plot_headshape(headshape,'vertexsize',5); view(0,0);
% plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
% axis off;
% grid on;
% set(gcf,'color','w')

%%
% figure
% ft_plot_vol(individual_headmodel, 'unit', 'mm');  %this is the brain shaped head model volume
% ft_plot_sens(t_data.all.grad, 'unit', 'mm', 'coilsize', 10);  %this is the sensor locations
% ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
% ft_plot_ortho(mri_realigned.anatomy, 'transform', mri_realigned.transform, 'style', 'intersect');

% quality check for the coregistration between headshape and mri
% headshape    = ft_convert_units(shape,    'mm');
% % headshapemri = ft_convert_units(shapemri, 'mm');
%
% v = headshapemri.pnt;
% f = headshapemri.tri;
% [f,v]=reducepatch(f,v, 0.2);
% headshapemri.pnt = v;
% headshapemri.tri = f;
%
% h = figure;
% subplot('position',[0.01 0.51 0.48 0.48]);hold on;
% % ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
% ft_plot_headshape(headshape,'vertexsize',5); view(180,-90);
% plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
% subplot('position',[0.51 0.51 0.48 0.48]);hold on;
% % ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
% ft_plot_headshape(headshape,'vertexsize',5); view(0,90);
% plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
% subplot('position',[0.01 0.01 0.48 0.48]);hold on;
% ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
% ft_plot_headshape(headshape,'vertexsize',5); view(90,0);
% plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
% subplot('position',[0.51 0.01 0.48 0.48]);hold on;
% % ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
% ft_plot_headshape(headshape,'vertexsize',5); view(0,0);
% plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
% axis off;
% grid on;
% set(gcf,'color','w')