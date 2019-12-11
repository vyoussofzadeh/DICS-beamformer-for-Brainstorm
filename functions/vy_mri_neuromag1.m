function [mri_realigned,individual_seg,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag1(cfg_main)

if exist(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']), 'file') == 2
    load(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']));
    load(fullfile(cfg_main.outputmridir,['mesh8mm_',cfg_main.subj,'.mat']));
    load(fullfile(cfg_main.outputmridir,['mesh10mm_',cfg_main.subj,'.mat']));
else
    if exist(cfg_main.mripfile,'file')== 2
        
        %%
        individual_mri = ft_read_mri(cfg_main.mripfile);
        %     ft_sourceplot([], individual_mri);
        individual_mri = ft_convert_units(individual_mri, 'mm');
        
        %%
        fid = cfg_main.fid;
        
        cfg = [];
        cfg.method = 'fiducial';
        cfg.coordsys = 'neuromag';
        cfg.fiducial.nas   = fid.SCS.NAS;
        cfg.fiducial.lpa   = fid.SCS.LPA;
        cfg.fiducial.rpa   = fid.SCS.RPA;
        cfg.spmversion     = 'spm12';
        mri_realigned = ft_volumerealign(cfg, individual_mri);
        
        %%
        cfg          = [];
        % cfg.method   = 'interactive';
        cfg.method = 'fiducial'; % the following voxel coords were determined interactive
        cfg.coordsys = 'spm';
        cfg.fiducial.ac   = fid.NCS.AC;
        cfg.fiducial.pc   = fid.NCS.PC;
        cfg.fiducial.xzpoint  = fid.NCS.IH;
        cfg.fiducial.right   = fid.SCS.RPA;
        mri_realigned     = ft_volumerealign(cfg, mri_realigned);
        
        %%
        %         ft_sourceplot([], mri_realigned); hold on
        %         plot3(fid.NAS(1,1), fid.NAS(1,2), fid.NAS(1,3), 'm*');
        %         plot3(fid.LPA(1,1), fid.LPA(1,2), fid.LPA(1,3), 'm*');
        %         plot3(fid.RPA(1,1), fid.RPA(1,2), fid.RPA(1,3), 'm*');
        
        %%
        headshape = ft_read_headshape(cfg_main.hsfile);
        headshape = ft_convert_units(headshape, 'mm');
        
        cfg = [];
        cfg.method = 'headshape';
        cfg.headshape.interactive = 'no';
        cfg.headshape.icp = 'yes';
        cfg.headshape.headshape = headshape;
        cfg.coordsys = 'neuromag';
        cfg.spmversion     = 'spm12';
        mri_realigned = ft_volumerealign(cfg, mri_realigned);
        
        
        %% == to check everything went right with co-reg!
        ft_determine_coordsys(mri_realigned, 'interactive', 'no')
        ft_plot_headshape(headshape);
        view([-90, 0]),
        
        %%
        cfg = [];
        cfg.output = 'brain';
        cfg.spmversion = 'spm12';
        individual_seg = ft_volumesegment(cfg, mri_realigned);
        individual_seg.transform = mri_realigned.transform;
        
        cfg = [];
        cfg.method = 'singleshell';
%         cfg.spmversion = 'spm12';
        individual_headmodel = ft_prepare_headmodel(cfg, individual_seg);
        
        %% Source modelling (warped with template)
        load temp_grid % low-res
        status = ft_hastoolbox('SPM12', 1, 0);
        cfg                 = [];
        cfg.grid.warpmni    = 'yes';
        cfg.spmversion     = 'SPM12';
        cfg.grid.nonlinear  = 'yes';
        cfg.grid.template   = template_grid;
        cfg.mri             = mri_realigned;
        cfg.grid.unit       = 'mm';
        individual_grid_10mm     = ft_prepare_sourcemodel(cfg);
        load temp_grid_8mm % high-res
        cfg.grid.template   = template_grid;
        individual_grid_8mm     = ft_prepare_sourcemodel(cfg);
        
        %%
        save(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']), 'mri_realigned','individual_seg','individual_headmodel','headshape');
        save(fullfile(cfg_main.outputmridir,['mesh8mm_',cfg_main.subj,'.mat']), 'individual_grid_8mm');
        save(fullfile(cfg_main.outputmridir,['mesh10mm_',cfg_main.subj,'.mat']), 'individual_grid_10mm');
    end
end

%% Quick inspection
if cfg_main.plotflag == 1
    figure;
    ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    hold on;
    ft_plot_headshape(headshape);
    ft_plot_mesh(individual_grid_10mm.pos(individual_grid_10mm.inside, :));
    view ([0 90])
end

%%
% create the individual grid from mni grid
%
%         temp = load('standard_sourcemodel3d8mm');
%         template_grid = ft_convert_units(template_grid, 'mm');
%         cfg = [];
%         cfg.grid.warpmni = 'yes';
%         %     cfg.grid.template = temp.sourcemodel;
%         cfg.grid.template = template_grid;
%         cfg.grid.nonlinear = 'yes';
%         cfg.mri = mri_realigned;
%         cfg.grid.unit     = 'mm';
%         individual_grid = ft_prepare_sourcemodel(cfg);
