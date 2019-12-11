function [sourcemodel] = vy_anatomy_sourcemodel2d_test(cfg, subject)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% if ischar(subject)
%   subject = streams_subjinfo(subject);
% end

anatomy_dir           = cfg.anatomy_dir;
inp_dir               = fullfile(anatomy_dir, subject);
sourcemodel_filename  =fullfile(anatomy_dir, [subject, '_sourcemodel.mat']); %string for saving the sourcemodel file


% load in the cortical sheet
filename = fullfile(inp_dir,['workbench/' subject, '.L.midthickness.8k_fs_LR.surf.gii']);
filename2 = strrep(filename, '.L.', '.R.');

sourcemodel = ft_read_headshape({filename, filename2});

% get the necessary coregistration information
datapath = fullfile(anatomy_dir);
load(fullfile(datapath,[subject,'_transform_vox']));
T1 = transform_vox;
load(fullfile(datapath,[subject,'_transform_vox2neuromag']));
T2 = transform_vox2neuromag;

sourcemodel = ft_transform_geometry((T2/T1), sourcemodel);
sourcemodel.inside = sourcemodel.atlasroi>0;
sourcemodel = rmfield(sourcemodel, 'atlasroi');

save(sourcemodel_filename, 'sourcemodel');

end

