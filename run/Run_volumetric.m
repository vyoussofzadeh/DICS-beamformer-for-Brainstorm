%% Volumetric-based analysis
mridir = fullfile(indir,subj,'brainstorm_db/anat');
d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
clear fid
if ~isempty(d)
    sMRI1 = d.name;
    load(sMRI1);
    fid.SCS = SCS;
    fid.NCS = NCS;
    mripfile = fullfile(mridir,'T1.nii');
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    cfg = [];
    cfg.megdata = t_data.pst.grad;
    cfg.mripfile = mripfile;
    cfg.hsfile = datafile; % headshape;
    cfg.fid = fid;
    cfg.outputmridir = outputmridir;
    cfg.subj = subj;
    cfg.plotflag = 2;
    cfg.atlas = atlas;
    cfg.indir = indir;
    cfg.outd.sub = outd.sub;
    cfg.flag = flag;
    outanat = vy_mri_neuromag2(cfg);
    %     vy_do_freesurfer(cfg);
end

%% Choosing mesh
if flag.warping == 1
    switch meshgridres
        case 1
            meshtag = 'lowres';
            %         load('standard_sourcemodel3d10mm');
            load temp_grid % low-res
            template_grid = ft_convert_units(template_grid, 'mm');
            individual_grid = outanat.individual_grid_10mm;
        case 2
            meshtag = 'highres';
            %         load('standard_sourcemodel3d8mm');
            load temp_grid_8mm % high-res
            individual_grid = outanat.individual_grid_8mm;
    end
else
    switch meshgrid
        case 1            
            individual_grid = outanat.individual_grid_10mm;
        case 2
            individual_grid = outanat.individual_grid_8mm;
    end
end

%% Anatomoy check!
saveflag = 2;
if flag.anatomycheck == 1
    
    cfg = [];
    cfg.saveflag = saveflag;
    cfg.headmodel = outanat.individual_headmodel;
    cfg.leadfield = individual_grid;
    cfg.mri_realigned  = outanat.mri_realigned;
    cfg.headshape = outanat.headshape;
    cfg.outputmridir = outputmridir;
    cfg.mtd = 'vol';
    vy_mri_inspection(cfg, t_data);
    %     vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
end

%%
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %

%%
switch method
    case 1
        %% LCMV
        %         mtag = 'lcmv_stat';
        mtag = 'lcmv';
        outd.vol = fullfile(outd.sub, mtag);
        cfg = [];
        cfg.allpath = allpath;
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        %         cfg.mtag = 'lcmv'; cfg.filterflag =  1;
        cfg.mtag = mtag; cfg.filterflag =  1;
        outd.vol = fullfile(outd.sub,cfg.mtag);
        %         cfg.fb   = [12,30]; % Hz
%         cfg.fb   = [16,25]; % Hz
        cfg.fb   = [14,27]; % Hz
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.toi       = {toi};
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        cfg.plotflag     = 1;
        switch mtag
            case 'lcmv'
                vy_source_lcmv(cfg, cln_data);
            case 'lcmv_stat'
                vy_source_lcmv_stats(cfg, cln_data);
        end
        
    case 2
        %%
        cfg = [];
        mtag = 'conn';
%         mtag = 'conn_bl'; cfg.fl = [1, 10];% band limited
%         mtag = 'conn_bs'; cfg.fb = 10; % band-stop
        outd.vol = fullfile(outd.sub,mtag);
        cfg.allpath = allpath;
        cfg.grid = individual_grid;
        cfg.headmodel = outanat.individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        cfg.mtag = mtag;
        cfg.toi  = toi;
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        switch mtag
            case 'conn'
                vy_network_light1(cfg,t_data) % conn-network analysis
            case {'conn_bs','conn_bl'}
                vy_network_light1(cfg, cln_data) % conn-network analysis
        end
        
    case 3
        %%
        mtag = 'dics';
%         mtag = 'dics_ratio';
%                 mtag = 'dics_stat';
        %         mtag = 'dics_fs';
        if rl == 1, mtag_lab = 'dics_res'; else, mtag_lab = 'dics'; end
        outd.vol = fullfile(outd.sub,mtag_lab);
        
        switch mtag
            
            case 'dics'
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
                cfg.freq_of_interest  = freq_of_interest; % Hz
                cfg.headmodel = outanat.individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                if flag.warping ==1
                cfg.template_grid = template_grid;
                end
                cfg.template_mri = template_mri;
                cfg.fmax = fmax;
                cfg.savedata = fullfile(outd.vol,[mtag,'_',subj]);
                cfg.flag = flag;
                vy_source_dics(cfg, ep_data);
                
            case 'dics_ratio'
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
                cfg.freq_of_interest   = freq_of_interest; % Hz
                cfg.headmodel = individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                cfg.template_grid = template_grid;
                cfg.template_mri = template_mri;
                cfg.savedata = fullfile(outd.vol,[mtag,'_',subj,'.mat']);
                vy_source_dics_ratio(cfg, ep_data);
                
            case 'dics_fs'
                %%
                anatomy_dir     = '/data/MEG/Clinical/ft_process/19/bednar_peggy/anat/BAK';
                load(fullfile(anatomy_dir,[subj,'_headmodel.mat']));
                load(fullfile(anatomy_dir,[subj,'_leadfield.mat']));
                load(fullfile(anatomy_dir,[subj,'_sourcemodel.mat']));
                
                if anatomy_check_flag == 1
                    cfg1 = [];
                    cfg1.saveflag = 2;
                    cfg1.headmodel = headmodel;
                    cfg1.sourcemodel = sourcemodel;
                    cfg1.leadfield = leadfield;
                    cfg1.mri = anatomy_dir;
                    cfg1.mtd = 'surf';
                    cfg1.headshape = headshape;
                    vy_mri_inspection(cfg1, t_data);
                    %                 individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag
                end
                cfg = [];
                cfg.headmodel = headmodel;
                cfg.sourcemodel = sourcemodel;
                cfg.leadfield   = leadfield;
                cfg.mtag = mtag;
                cfg.sens = sens;
                cfg.subj = subj;
                cfg.outputdir = outd.vol;
                vy_source_dics_fs(cfg, ep_data);
                
            case 'dics_stat'
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
                cfg.f   = 20; % Hz
                cfg.headmodel = individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                cfg.template_grid = template_grid;
                cfg.template_mri = template_mri;
                vy_source_dics_stats(cfg, ep_data);
        end
        
end
