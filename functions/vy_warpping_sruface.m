function vy_warpping_surface()


%-
% template_mri = ft_read_mri('F:\My Matlab\MEG\HCP\megconnectome-3.0\template\T1.nii'); % based on mni151
% template_mri.coordsys = 'spm';  % inform fieldtrip that 'spm' coordsystem
% 
% cfg = [];
% cfg.spmversion = 'spm12';
% template_seg = ft_volumesegment(cfg, template_mri); % segment the template mri
% save('seg','template_seg');
load seg

cfg = [];
cfg.method = 'singleshell';
template_headmodel = ft_prepare_headmodel(cfg, template_seg);   % create headmodel from segmented template
template_headmodel = ft_convert_units(template_headmodel, 'mm'); % convert to cm (CTF native coordinate system)

cfg = [];
% cfg.grid.resolution = 8;
% cfg.grid.tight = 'yes';
% cfg.inwardshift = -8;                        % include region just outside head as 'inside', to avoid rim of nonactivation
cfg.headmodel = template_headmodel;
template_grid = ft_prepare_sourcemodel(cfg);        % create the template grid (coords of interest)

figure
hold on
ft_plot_vol(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% ft_plot_mesh(template_grid.pos(template_grid1.inside,:));
ft_plot_mesh(template_grid.pos);


% save('temp_grid_8mm','template_grid');
save('temp_grid_surface','template_grid');

% load temp_grid_8mm

%%
% template_grid = ft_convert_units(template_grid, 'mm');
% cfg                 = [];
% cfg.grid.warpmni    = 'yes';
% cfg.grid.nonlinear  = 'yes';
% % cfg.grid.templatemri = template_mri;
% cfg.grid.template   = template_grid;
% cfg.mri        = mri_realigned;
% cfg.grid.unit     = 'mm';
% individual_grid    = ft_prepare_sourcemodel(cfg);
% 
% %%
% % individual_headmodel.coordsys = 'bti';
% % individual_headmodel  = ft_convert_units(individual_headmodel);
% 
% % sourcemodel3d = ft_convert_units(sourcemodel3d, 'mm'); % convert to cm (CTF native coordinate system)
% figure;
% % subplot(1, 2, 1);
% ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');
% alpha 0.5;
% camlight;
% hold on;
% ft_plot_mesh(individual_grid.pos(individual_grid.inside,:));

