function vy_network_light1(cfg_main, data_in)

%% BPF

%%
mtag = cfg_main.mtag;
switch mtag
    case 'conn_bl' % conn band limited
        cfg = [];
        cfg.hpfilter = 'yes';
        cfg.lpfilter = 'yes';
        cfg.hpfreq = cfg_main.fb(1);
        cfg.lpfreq = cfg_main.fb(2);
        fcln_data = ft_preprocessing(cfg, data_in);
        hp = cfg_main.fb(2);
        lp = cfg_main.fb(2);
        %-
        cfg = [];
        cfg.savefile = []; cfg.saveflag = 2; cfg.foilim = [2 40];
        cfg.plotflag  = 1; cfg.tapsmofrq = 5; cfg.taper     = 'hanning';
        vy_fft(cfg, fcln_data);
        %-
        fep_data = vy_epoch(fcln_data, cfg_main.toi);
        cfg = [];
        fep_data.app = ft_appenddata(cfg,fep_data.bsl,fep_data.pst);
        data_in = vy_timelock(fep_data);

    case 'conn_bs' % conn band limited
        
        cfg = [];
        cfg.bsfilter = 'yes';
        %     cfg.bsfreq = [29 32]; % or whatever you deem appropriate
        cfg.bsfreq = [cfg_main.fb-2 cfg_main.fb+2]; % or whatever you deem appropriate
        %     cfg.bsfreq = [8 12;29 32]; % or whatever you deem appropriate
        fcln_data = ft_preprocessing(cfg, data_in);
        
        cfg = [];
        cfg.savefile = [];
        cfg.saveflag = 2;
        cfg.foilim = [2 40];
        cfg.plotflag  = 1;
        cfg.tapsmofrq = 5;
        cfg.taper     = 'hanning';
        vy_fft(cfg, fcln_data);
        grid on
        grid minor
        title(['After band-stop filtering-',cfg_main.subj]);
        %         cfg.hpfilter = 'yes';
        %         cfg.lpfilter = 'yes';
        %         cfg.hpfreq = cfg_main.fb(1);
        %         cfg.lpfreq = cfg_main.fb(2);
        %         fcln_data = ft_preprocessing(cfg, data_in);
        %         hp = cfg_main.fb(2);
        %         lp = cfg_main.fb(2);
        %-
        fep_data = vy_epoch(fcln_data, cfg_main.toi);
        cfg = [];
        fep_data.app = ft_appenddata(cfg,fep_data.bsl,fep_data.pst);
        data_in = vy_timelock(fep_data);
end

%%

cfg = [];
cfg.grid = cfg_main.grid;
cfg.headmodel = cfg_main.headmodel;
cfg.mtd = 'lcmv';
% cfg.mtd = 'sam';
s_data = vy_source(cfg, data_in);
% s_data = vy_source(t_data, cfg_main.grid, cfg_main.headmodel);

%%
% using older ver of ft for network analysis
% restoredefaultpath
% addpath(genpath(cfg_main.allpath.ft_oldp.ft_old));
% addpath(genpath(cfg_main.p.hcp_path));
% addpath(genpath(cfg_main.p.cd_org));
%
restoredefaultpath
addpath(genpath(cfg_main.allpath.ft_old));
addpath(genpath(cfg_main.allpath.hcp_path));
addpath(genpath(cfg_main.allpath.cd_org));

%%
mtd = 'plv';
mtd_par = 'plvspctrm';
gtm = 'eigenvector_cent';

conn_par = [];
conn_par.method   = mtd;
conn_par.idx      = mtd_par;
conn_par.complex  = [];

net_par.gtm       = gtm ; in2 = 2;
net_par.threshold = 0.7;
[source_conn_bsl, network_bsl] = vy_conn(s_data.bsl,conn_par,net_par);
[source_conn_pst, network_pst] = vy_conn(s_data.pst,conn_par,net_par);

%%
cfg = [];
cfg.parameter = gtm;
cfg.operation = 'x1-x2';
network_diff_lcmv = ft_math(cfg,network_pst,network_bsl);
network_diff_lcmv.pos     = cfg_main.template_grid.pos;
network_diff_lcmv.dim     = cfg_main.template_grid.dim;
network_diff_lcmv.inside  = cfg_main.template_grid.inside;
network_diff_lcmv.(gtm)   = zscore(network_diff_lcmv.(gtm));
network_diff_lcmv.eigenvector_cent(network_diff_lcmv.eigenvector_cent<0)=0;

%% Revert to new ft!
restoredefaultpath
addpath((cfg_main.allpath.ft_path));
ft_defaults
addpath(genpath(cfg_main.allpath.hcp_path));
addpath(genpath(cfg_main.allpath.cd_org));
addpath(genpath(cfg_main.allpath.exfig_path));

%%
toi = cfg_main.toi;
if size(toi,1) > 1
    sname2 = [num2str(toi(2,1)),'_',num2str(toi(2,2)),'sec'];
    disp([num2str(toi(2,1)),'_',num2str(toi(2,2)),'sec is analysing'])
else
    sname2 = [num2str(toi(1)),'_',num2str(toi(2)),'sec'];
    disp([num2str(toi(1)),'_',num2str(toi(2)),'sec is analysing'])
end

%%

savepath = fullfile(cfg_main.outputdir);
if exist(savepath, 'file') == 0, mkdir(savepath), end
savedata = fullfile(savepath,['n_',sname2,'_',cfg_main.subj,'.mat']);
save(savedata, 'network_diff_lcmv', '-v7.3');

mtd = 'network_evc';

% oldmask = network_diff_lcmv.(gtm);
% network_diff_lcmv.mask = (oldmask - min(oldmask(:))) ./ (max(oldmask(:)) - min(oldmask(:))); %
% % set the new mask to range between 0 and 1
% network_diff_lcmv.mask(oldmask < 0.9) = 0; % mask out non-significant voxels
% savefig = fullfile(outputdir_dics,[num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);
% savefig = fullfile(savepath,[mtd,'_1_',cfg_main.subj]);
savefig = fullfile(savepath,['n_',sname2,'_',cfg_main.subj]);

cfg = [];
cfg.mask = gtm;
cfg.loc = 'max';
cfg.template = cfg_main.template_mri;
cfg.savefile = savefig;
cfg.volnorm     = 2; % yes: 1
network_int = vy_source_plot(cfg, network_diff_lcmv);

% vy_mapvisualisation(network_int_lcmv,gtm,0.6,savep)
% vy_mapvisualisation(network_int,gtm,0.6, []);


clear savepath
savepath{1} = fullfile(cfg_main.outputdir,['_R_',cfg_main.subj]);
savepath{2} = fullfile(cfg_main.outputdir,['_L_',cfg_main.subj]);
% vy_mapvisualisation(network_int,gtm,0.6, savepath);
% vy_mapvisualisation(network_int,gtm, 0.6, []);

cfg = [];
cfg.subj = cfg_main.subj;
cfg.mask = gtm;
cfg.thre = 0.6;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, network_int);

% cfg = [];
% cfg.subj = cfg_main.subj;
% cfg.mask = gtm;
% cfg.thre = 0.6;
% cfg.savepath = savepath;
% cfg.savepath = [];
% vy_mapvisualisation(cfg, network_int);

% savenii = fullfile(savepath,['n_',subj,'.nii']);
% vy_savenifti(network_int_lcmv, gtm, savenii);

%%


% param = [];
% param.mask = gtm;
% param.loc = 'max';
% network_int_lcmv = vy_source_plot(network_diff_lcmv,cfg_main.template_mri,param,2);
% hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
% saveformat = '-png';
% pixdim     = '-m8';
% export_fig(savefig, saveformat, pixdim)
% print(gcf, '-dpdf',[savefig,'.pdf']);

% clear savep
% savep{1} = fullfile(savepath,[mtd,'_2_',cfg_main.subj]);
% savep{2} = fullfile(savepath,[mtd,'_3_',cfg_main.subj]);

% cfg = [];
% cfg.maskparam = gtm;
% cfg.save.savepath =  savep;
% % cfg.saveformat = '-eps';
% cfg.save.saveformat = '-png';
% cfg.save.pixdim     = 12;
% cfg.projthresh      = 0.6;
% vy_surfmap(cfg, network_int_lcmv);
% cfg = [];
% cfg.funparameter = param.mask;
% cfg.method = 'ortho';  % plot slices
% ft_sourceplot(cfg, network_diff_lcmv);


%%
% network_int1 = network_int_lcmv;
% thre = 0.5;
% network_int1.eigenvector_cent(network_int1.eigenvector_cent < thre*max(network_int1.eigenvector_cent(:)))=0;
% % network_int1.eigenvector_cent(network_int1.eigenvector_cent > 0) = 1;
% cfg = [];
% cfg.filetype  = 'nifti';
% cfg.parameter = gtm;
% % cfg.filename  = './output';
% cfg.filename  = fullfile(outputdir_net,'surf.nii');
% ft_volumewrite(cfg, network_int1)
%
% restoredefaultpath
% addpath(genpath([cd_org,'/functions']));
% % close all
% addpath(connpath);
% addpath(genpath(spm_path))
% h = get(0, 'Children');
% if isempty(findobj(h,'tag','CONN functional connectivity toolbox'))
%     conn
% end
%
% filenameVOL = [];
% FSfolder = fullfile(connpath,'/utils/surf');
% sphplots = [];
% connplots = [];
% facealpha = 1;
% position = [-1 0  0];
% conn_mesh_display(fullfile(outputdir_net,'surf.nii'),[],FSfolder);

%% revert to the newer ft!
% restoredefaultpath
% addpath(genpath(cfg_main.p.ft_path));
% addpath(genpath(cfg_main.p.hcp_path));
% addpath(genpath([cfg_main.p.cd_org,'/functions']));
% addpath(genpath([cfg_main.p.cd_org,'/Data_file']));

%% parcellation - aal (132 rois)

% network_diff_lcmv.eigenvector_cent = (network_diff_lcmv.eigenvector_cent)./max(network_diff_lcmv.eigenvector_cent);
% atlas = ft_read_atlas('/data/MEG/Vahab/Github/MCW-MEGlab/tools/ft_packages/fieldtrip_20190419/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');
% [~, data_intpar, coor] = vy_parcellate(network_diff_lcmv, cfg_main.atlas, gtm); % this part works with ft 2018 only!
% data_intpar.eigenvector_centdimord = 'chan';
%

%% Conn
% conn_par.conn_thre = 0.95;
% conn_ratio = vy_connvis(source_conn_pst,source_conn_bsl,conn_par, individual_headmodel, network_diff_lcmv);
% view([156,47]);

%
%% ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, gtm);
% disp(ROI_sel)
% savepath = fullfile(outputdir_net,['n_ROIs_',subj]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);


end