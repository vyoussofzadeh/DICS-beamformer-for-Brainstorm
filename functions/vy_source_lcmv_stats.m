function vy_source_lcmv_stats(cfg_main, data_in)

outputdir_lcmv = cfg_main.outputdir;
if exist(outputdir_lcmv, 'file') == 0, mkdir(outputdir_lcmv), end

sname3 = cfg_main.subj;

mtag = cfg_main.mtag;
fflag = cfg_main.filterflag;

for i=1:length(cfg_main.toi)
    
    toi = cfg_main.toi{i};
    sname2 = [num2str(toi(2,1)),'_',num2str(toi(2,2)),'sec'];
    
    disp('============');
    disp([num2str(toi(2,1)),'_',num2str(toi(2,2)),'sec is analysing'])
    disp('============');
    
    if fflag == 1
        
        cfg = [];
        cfg.hpfilter = 'yes';
        cfg.lpfilter = 'yes';
        cfg.hpfreq = cfg_main.fb(1);
        cfg.lpfreq = cfg_main.fb(2);
        fcln_data = ft_preprocessing(cfg, data_in);
        hp = cfg.hpfreq;
        lp = cfg.lpfreq;
        sname1 = [num2str(hp),'_',num2str(lp),'Hz'];
        
        %%
        fep_data = vy_epoch(fcln_data, toi);
        cfg = [];
        fep_data.app = ft_appenddata(cfg,fep_data.bsl,fep_data.pst);
        t_data = vy_timelock(fep_data);
        savefig = fullfile(outputdir_lcmv,[sname1,'_',sname2,'_1_',sname3]);
        
    else
        t_data = cfg_main.t_data;
        savefig = fullfile(outputdir_lcmv,[sname2,'_1_',sname3]);
    end
    
    %%
    cfg = [];
    cfg.grid = cfg_main.grid;
    cfg.headmodel = cfg_main.headmodel;
    cfg.mtd = mtag;
    s_data_lcmv = vy_source(cfg, t_data);
    
    %%
    %     if cfg_main.plotflag == 1
    
    stat = vy_source_stat_montcarlo(s_data_lcmv);
    stat.pos     = cfg_main.template_grid.pos;
    stat.dim     = cfg_main.template_grid.dim;
    stat.inside  = cfg_main.template_grid.inside;
    
    tmp = stat.stat;
    tmp2 = zeros(size(stat.pos,1),1);
    tmp2(stat.inside) = tmp;
    
    stats1  = stat;
    stats1.stat =  tmp2;
    stats1.mask = stat.inside;
    stats1.stat(stats1.stat>0)=0;
    stats1.stat(isnan(stats1.stat))=0;
    
    %         savefig = fullfile(outputdir_lcmv,[num2str(hp),'_',num2str(lp),'Hz','_stats_1_',cfg_main.subj]);
    if cfg_main.plotflag == 1
        
        cfg = [];
        cfg.mask = 'stat';
        cfg.loc  = 'min';
        cfg.template = cfg_main.template_mri;
        cfg.savefile = savefig;
        cfg.volnorm     = 2; % yes: 1
        source_int_lcmv = vy_source_plot(cfg, stats1);
        
        vy_mapvisualisation(source_int_lcmv,cfg.mask,0.6, []);
    end
    
    stats{i} =  stats1.stat;
    
end

save([cfg_main.outputdir, '/', [cfg_main.mtag,'_',sname3,'.mat']], 'stats1', 'stats','-v7.3');
% plot avg
mm  = mean(cell2mat(stats),2);
stats1.stat = mm;
savefig = fullfile(outputdir_lcmv,['avg_slice_',sname3]);

cfg = [];
cfg.mask = 'stat';
cfg.loc = 'min';
cfg.template = cfg_main.template_mri;
cfg.savefile = savefig;
cfg.volnorm     = 2; % yes: 1
source_int_lcmv = vy_source_plot(cfg, stats1);

clear savepath
savepath{1} = fullfile(outputdir_lcmv,['avg_R_',sname3]);
savepath{2} = fullfile(outputdir_lcmv,['avg_L_',sname3]);
vy_mapvisualisation(source_int_lcmv,cfg.mask,0.6, savepath);
%         vy_mapvisualisation(source_int_lcmv,cfg.mask,0.6, []);

end
%%
% save nii
% savenii = fullfile(outputdir_lcmv,['s_lcmv_',subj,'_',run,'.nii']);
% vy_savenifti(source_int_lcmv,'pow',savenii);
% vy_savenifti(source_int_lcmv,'pow','test.nii');

%% parcellation
% [~, data_intpar, coor] = vy_parcellate(source_diff_lcmv, atlas,'pow');
%
% savepath = fullfile(outputdir_lcmv,[mtd,'_par_1',subj,'_',run,'.mat']);
% save(savepath, 'data_intpar', '-v7.3');
%
% clear savepath
% savepath{1} = fullfile(outputdir_lcmv,[mtd,'_par_2',subj,'_',run]);
% savepath{2} = fullfile(outputdir_lcmv,[mtd,'_par_3',subj,'_',run]);
% vy_mapvisualisation(data_intpar,'pow',0.6, savepath)
%
% % save nii
% savenii = fullfile(outputdir_lcmv,['s_lcmv_par',subj,'_',run,'.nii']);
% vy_savenifti(data_intpar,'pow',savenii);
%
% % ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, 'pow');
% disp(ROI_sel)
% savepath = fullfile(outputdir_lcmv,['s_dics_ROIs_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);