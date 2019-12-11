addpath(allpath.connpath);
addpath(allpath.spm_path);

%%
group_source = ft_read_mri([tsk,'_groupave.nii']);

% Opt = [];
% Opt.savenii = 0; Opt.savefig = 0;
% Opt.savename = [tsk,'_groupave1'];
% vy_surfce_vis2(group_source,[tsk,'_groupave.nii'], Opt);

%%
projthresh = 0.80;
s_vol = vy_vol_thresh(group_source, projthresh, 'anatomy'); % abs

Opt = [];
Opt.savenii = 1; Opt.savefig = 0;
Opt.savename = [tsk,'_groupave_thre'];
vy_surfce_vis2(s_vol,[tsk,'_groupave_thre.nii'], Opt);

%%
% view([110,22])
% view([-110,22])

