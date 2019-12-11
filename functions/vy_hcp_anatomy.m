function [individual_grid, mri_realigned, individual_seg] = vy_hcp_anatomy(dicomfile,hsfile, outputdir,subjectid)

% load temp_grid

% sfile = fullfile(outputdir,[subjectid,'_MEG_anatomy_realigned.mat']);
% 
% if exist(sfile, 'file') == 2
%     load(sfile);
%     load(fullfile(outputdir,[subjectid,'_MEG_anatomy_headshapemri.mat']));
%     load(fullfile(outputdir,[subjectid,'_MEG_anatomy_headshape.mat']));
% else
% structuralpreprocdir = [];
dopipeinteractive = 1;
dopipeautomatic   = 1;
doqualitycheck    = 1;

dosourcemodel2d = 0;
    
    hcp_anatomy
% end




%
%     individual_mri = ft_read_mri(mripath);
%     individual_mri.coordsys = '4d';
%     ft_sourceplot([], individual_mri);
%
%     cfg = [];
%     cfg.method = 'flip';
%     individual_mri = ft_volumereslice(cfg, individual_mri);
%     individual_mri = ft_convert_units(individual_mri, 'm');
%     ft_sourceplot([], individual_mri);
%
%
%     headshape = ft_read_headshape('F:\UTHSC\VNS\MEG\C-101\CRM\1\hs_file');
%     mri_fids = headshape.fid.pos;
%
%     % the MRI is neither expressed in MNI, nor in Neuromag coordinates
%     ft_determine_coordsys(individual_mri, 'interactive', 'no');
%     hold on; % add the subsequent objects to the same figure
%     %     ft_plot_headshape(headshape);
%
%     plot3(mri_fids(1,1), mri_fids(1,2), mri_fids(1,3), 'm*');
%     plot3(mri_fids(2,1), mri_fids(2,2), mri_fids(2,3), 'm*');
%     plot3(mri_fids(3,1), mri_fids(1,2), mri_fids(3,3), 'm*');
%
%
%     %     cfg = [];
%     %     cfg.fiducial.nas = vox_fids(1,:);
%     %     cfg.fiducial.lpa = vox_fids(2,:);
%     %     cfg.fiducial.rpa = vox_fids(3,:);
%     %     cfg.coordsys = '4d';
%     %     mri_realigned = ft_volumerealign(cfg, individual_mri);
%
%
%     cfg = [];
%     cfg.method = 'interactive';
%     cfg.coordsys = '4d';
%     mri_realigned = ft_volumerealign(cfg, individual_mri);
%
%
% %     cfg = [];
% %     cfg.method = 'headshape';
% %     cfg.coordsys = '4d';
% %     cfg.headshape = headshape;
% %     mri_realigned2 = ft_volumerealign(cfg, mri_realigned);
%
%     % create the individual grid from mni grid
%     cfg = [];
%     cfg.grid.warpmni = 'yes';
%     cfg.grid.template = template_grid;
%     cfg.grid.nonlinear = 'yes';
%     cfg.mri = mri_realigned;
%     individual_grid = ft_prepare_sourcemodel(cfg);
%     %

%
%     save(fullfile(p,['ana_',subj,'.mat']), 'mri_realigned','individual_seg','individual_grid');
%
% end
