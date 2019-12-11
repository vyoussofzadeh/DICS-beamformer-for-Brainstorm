
outputdir1 = fullfile(outd.sub, 'spm_source');
if exist(outputdir1, 'file') == 0
    mkdir(outputdir1);   % create a directory
end
addpath(genpath(allpath.spm_path))
cfg = [];
cfg.toilim = [-0.4 1.2];
eint_data = ft_redefinetrial(cfg, cln_data);

%%
if exist(mripfile, 'file') == 2
    cd(outputdir1);
    cfg = [];
    cfg.allpath = allpath;
    cfg.datafile = datafile;
    cfg.eint_data = eint_data;
    cfg.sub       = subj;
    cfg.mripfile = mripfile;
    vy_forward_spm_meg(cfg);
    cd(outputdir1)
end
cd ..

% revert to the newer ft!
restoredefaultpath
addpath((allpath.ft_path));
ft_defaults
addpath(genpath(allpath.hcp_path));
addpath(genpath(allpath.cd_org));