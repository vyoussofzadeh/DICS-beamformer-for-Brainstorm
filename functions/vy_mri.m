function [mri_realigned,individual_seg,individual_headmodel,headshape] = vy_mri(mripath,hsfile,outputmridir,subj)

if exist(fullfile(outputmridir,['ana_',subj,'.mat']), 'file') == 2
    load(fullfile(outputmridir,['ana_',subj,'.mat']),'mri_realigned','individual_seg','individual_grid','individual_headmodel','headshape');
else   
    if exist(mripath,'file')== 2
        
        individual_mri = ft_read_mri(mripath);
        individual_mri.coordsys = 'spm';
        %     figure, ft_sourceplot([], individual_mri);
        individual_mri = ft_convert_units(individual_mri, 'mm');
        
        headshape = ft_read_headshape(hsfile);
        headshape = ft_convert_units(headshape, 'mm');
        
        fprintf('Please identify the LPA, RPA, nasion, and a point on the positive Z-axis\n');
        cfg = [];
        cfg.method = 'interactive';
        cfg.coordsys = 'neuromag';
        individual_mri = ft_volumerealign(cfg, individual_mri);
        
        cfg = [];
        cfg.method = 'headshape';
        cfg.headshape.interactive = 'yes';
        cfg.headshape.icp = 'yes';
        cfg.headshape.headshape = headshape;
        cfg.coordsys = 'neuromag';
        mri_realigned = ft_volumerealign(cfg, individual_mri);
        if  det(mri_realigned.transform(1:3,1:3))>0
            disp('transformation-matrix is right-handed: right is right')
        else
            disp('transformation-matrix is left-handed: right is left')
        end
        
        ft_determine_coordsys(mri_realigned, 'interactive', 'no')
        ft_plot_headshape(headshape);
        pause % to check everything went right with co-reg!
        
        cfg = [];
        cfg.output = 'brain';
        cfg.spmversion = 'spm12';
        individual_seg = ft_volumesegment(cfg, mri_realigned);
        individual_seg.transform = mri_realigned.transform;       
        
        cfg = [];
        cfg.method = 'singleshell';
        individual_headmodel = ft_prepare_headmodel(cfg, individual_seg);
        
        % create the individual grid from mni grid
        
        %     temp = load('F:\My Matlab\MEG\HCP\megconnectome-3.0\template\standard_sourcemodel3d8mm');
%         template_grid = ft_convert_units(template_grid, 'mm');
%         cfg = [];
%         cfg.grid.warpmni = 'yes';
%         %     cfg.grid.template = temp.sourcemodel;
%         cfg.grid.template = template_grid;
%         cfg.grid.nonlinear = 'yes';
%         cfg.mri = mri_realigned;
%         cfg.grid.unit     = 'mm';
%         individual_grid = ft_prepare_sourcemodel(cfg);
        
        
        save(fullfile(outputmridir,['ana_',subj,'.mat']), 'mri_realigned','individual_seg','individual_headmodel','headshape');
    end
end