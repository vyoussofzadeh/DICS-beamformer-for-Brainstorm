%         if exist(savepath_sourcemodel, 'file') == 2
%             load(savepath_sourcemodel)
%         else
addpath(genpath(spm_path))
d = ['spm_source\m',subj];
D = spm_eeg_load(d);
datareg = D.inv{1}.datareg;
mesh = spm_eeg_inv_transform_mesh(datareg.fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
allmeshvert_mni = D.inv{1}.mesh.tess_mni.vert;

hsfile = fullfile(datafolder1,'hs_file'); % headshape
headshape = ft_read_headshape(hsfile);
headshape = ft_convert_units(headshape, 'mm');

% create the source model
sourcemodel = [];
sourcemodel.pos = mesh.tess_ctx.vert;
sourcemodel.tri = mesh.tess_ctx.face;
sourcemodel.inside = true(size(sourcemodel.pos,1),1);
sourcemodel.unit = 'mm';

% create the head model
hdm                = [];
hdm.bnd.pos        = mesh.tess_iskull.vert; % template's inner skull surface
hdm.bnd.tri        = mesh.tess_iskull.face;
hdm.type           = 'singleshell';
hdm.unit           = 'mm';

% compute the leadfield
cfg             = [];
cfg.senstype    = 'MEG';
cfg.grid        = sourcemodel;
cfg.headmodel   = hdm;
cfg.channel     = {'MEG'};
lf              = ft_prepare_leadfield(cfg, t_data.app);

individual_grid = lf;
individual_headmodel = hdm;

%     figure;
%     ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
%     hold on;
%     ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
%     view ([-10 40 10])

% save(savepath_sourcemodel, 'individual_headmodel','individual_grid', 'allmeshvert_mni','-v7.3');
%         end
%- anatomy inspections ( headmodel, mesh, alignment, ...)
%  vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,sens, outputmridir)
