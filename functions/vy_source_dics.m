function vy_source_dics(cfg_main, ep_data)

if exist(cfg_main.savedata, 'file') == 2
    disp('source has already been computed, skipping')
    disp('Yes = 1')
    disp('No = 2');
    %     ask = input('calculating again?');
    ask = 1;
else
    ask = 1;
end
if ask==1
    
    Run_fft_4dics
    %%
    cfg = [];
    cfg.headmodel = cfg_main.headmodel;
    % cfg.sourcemodel = cfg_main.sourcemodel;
    cfg.grid = cfg_main.grid;
    cfg.mtag = cfg_main.mtag;
    s_data_dics = vy_source_freq(cfg, f_data);
    % s_data_dics = vy_source_freq(f_data, cfg_main.grid, cfg_main.headmodel, cfg_main.mtag);
    
    %%
    cfg = [];
    cfg.parameter = 'pow';
    cfg.operation = '(x1-x2)/(x1+x2)';
    source_diff_dics = ft_math(cfg,s_data_dics.pst,s_data_dics.bsl);
    source_diff_dics.pos     = cfg_main.template_grid.pos;
    source_diff_dics.dim     = cfg_main.template_grid.dim;
    source_diff_dics.inside  = cfg_main.template_grid.inside;
    source_diff_dics.pow(source_diff_dics.pow>0)=0;
    
    %%
    % outputdir_dics = fullfile(outputdir,'dics');
    % if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end
    toi = cfg_main.toi;
    save([cfg_main.savedata,'_',num2str(f),'Hz.mat'], 'source_diff_dics', 'f','toi','-v7.3');
    
    savefig = fullfile(outputdir_dics,[num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);    
    if size(toi,1) == 1
        savefig = fullfile(outputdir_dics,[num2str(toi(1)),'_',num2str(toi(2)),'sec_',num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);
    end
    
    cfg = [];
    cfg.mask = 'pow';
    cfg.loc = 'min';
    cfg.template = cfg_main.template_mri;
    cfg.savefile = savefig;
    cfg.volnorm     = 2; % yes: 1
    source_int_dics = vy_source_plot(cfg, source_diff_dics);
    set(gcf,'name',cfg_main.subj,'numbertitle','off')
        
    clear savepath
    savepath{1} = fullfile(outputdir_dics,[num2str(f),'Hz','_2_',cfg_main.subj]);
    savepath{2} = fullfile(outputdir_dics,[num2str(f),'Hz','_3_',cfg_main.subj]);
    
    cfg = [];
    cfg.subj = cfg_main.subj;
    cfg.mask = 'pow';
    cfg.thre = 0.6;
    cfg.savepath = savepath;
    vy_mapvisualisation(cfg, source_int_dics);
    % vy_mapvisualisation(source_int_dics,cfg.mask,0.6, []);
    
else
    disp('source estimatetion was not estimated!')
end
%%
% mtd = 'source_dics';
% param = [];
% param.mask = 'pow';
% param.loc = 'min';
% source_int_dics = vy_source_plot(source_diff_dics,cfg_main.template_mri,param,2);
% % savefig = fullfile(outputdir_dics,[num2str(f),'Hz','_1_',cfg_main.subj]);
% % hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
%
% clear savepath
% savepath{1} = fullfile(outputdir_dics,[num2str(f),'Hz','_2_',cfg_main.subj]);
% savepath{2} = fullfile(outputdir_dics,[num2str(f),'Hz','_3_',cfg_main.subj]);
% vy_mapvisualisation(source_int_dics,'pow',0.6, savepath);
% % vy_mapvisualisation(source_int_dics,'pow',0.6, []);

% restoredefaultpath
% addpath((cfg_main.allpath.ft_path));
% ft_defaults
% addpath(genpath(cfg_main.allpath.hcp_path));
% addpath(genpath(cfg_main.allpath.cd_org));
% addpath(genpath(cfg_main.allpath.exfig_path));
%
% cfg = [];
% cfg.maskparam = 'pow';
% cfg.save.savepath =  savepath;
% % cfg.saveformat = '-eps';
% cfg.save.saveformat = '-png';
% cfg.save.pixdim     = 12;
% cfg.projthresh      = 0.6;
% vy_surfmap(cfg, source_int_dics);

% save nii
% savenii = fullfile(outputdir_dics,['s_dics_',subj,'_',run,'.nii']);
% vy_savenifti(source_int_dics,'pow',savenii);

%% parcellation
% [~, data_intpar, coor] = vy_parcellate(source_diff_dics, atlas,'pow');
%
% savepath = fullfile(outputdir_dics,[mtd,'_par_1',subj,'_',run,'.mat']);
% save(savepath, 'data_intpar', '-v7.3');
%
% clear savepath
% savepath{1} = fullfile(outputdir_dics,[mtd,'_par_2',subj,'_',run]);
% savepath{2} = fullfile(outputdir_dics,[mtd,'_par_3',subj,'_',run]);
% vy_mapvisualisation(data_intpar,'pow',0.6, savepath)
%
% % save nii
% savenii = fullfile(outputdir_dics,['s_dics_par',subj,'_',run,'.nii']);
% vy_savenifti(data_intpar,'pow',savenii);
%
% % ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, 'pow');
% disp(ROI_sel)
% savepath = fullfile(outputdir_dics,[mtd,'_par_roi',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

%%
% load('cortex_inflated_shifted.mat')
% % ft_sourceplot configuration
% cfg               = [];
% cfg.funparameter  = 'pow';
% cfg.maskparameter = 'pow';
% cfg.maskstyle     = 'colormix';
% cfg.method        = 'surface';
% cfg.funcolormap   = flipud(brewermap(65, 'RdBu'));
% cfg.camlight      = 'no';
% cfg.colorbar      = 'yes';
%
% opacityratio      = 0.5;
%
% % general fig
% cfg.funcolorlim   = 'maxabs';
%
% sp        = keepfields(source_int_dics, {'pow' 'time' 'dimord'});
% % sp.stat   = mean(s.stat.*double(s.posclusterslabelmat==1), 3);
% sp.dimord = 'pos_time';
% sp.time   = 0;
% sp.pos    = ctx.pos;
% sp.tri    = ctx.tri;
% % Time slices
% % times           = [1 2 3 4];
% maxabs          = max(abs(source_int_dics.pow(:)));
% cfg.funcolorlim = [-maxabs maxabs];
%
% % lateral side
% cfg.opacitylim = cfg.funcolorlim.*opacityratio;
% ft_sourceplot(cfg, source_int_dics);
%
% view([90 0]); camlight; material dull; title('test');

%%
% m.tri = ctx.tri;
% m.pos = ctx.pos;
%
%
% figure;
% set(gcf,'color','w');
% ft_plot_mesh(m, 'vertexcolor', double(ba42)); colormap parula
% view([110 10]);
% l = camlight;
% material dull;

%%
% cfg              = [];
% cfg.method       = 'mtmfft';
% cfg.output       = 'pow';
% cfg.pad          = 'nextpow2';
% cfg.output       = 'fourier';
% cfg.keeptrials   = 'yes';
% cfg.foilim       = [f,f];
% cfg.tapsmofrq    = 4;
% cfg.taper        = 'hanning';
% % cfg.taper        = 'dpss';
% cfg.pad          = 4;
% cfg.keeptrials   = 'yes';
% freq_pst         = ft_freqanalysis(cfg, ep_data.pst);
%
%
% sourceDiff = s_data_dics.pst;
% sourceDiff.avg.pow    = (s_data_dics.pst.avg.pow - s_data_dics.bsl.avg.pow) ./ s_data_dics.bsl.avg.pow;
%
% [maxval, maxpowindx] = max(sourceDiff.avg.pow);
% sourceDiff.pos(maxpowindx, :)
%
% cfg = [];
% cfg.metheod = 'dics';
% cfg.dics.lambda = '5%';
% cfg.frequency    = f_data.app.freq;
% cfg.grid = cfg_main.grid;
% cfg.headmodel = cfg_main.headmodel;
% cfg.dics.keepfilter = 'yes';
% cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
% sourceavg = ft_sourceanalysis(cfg, f_data.app);
% cfg.grid.filter = sourceavg.avg.filter;
%
% virtualgrid = rmfield(cfg.grid, 'filter');
% virtualgrid.pos = virtualgrid.pos(maxpowindx,:);
% virtualgrid.inside = virtualgrid.inside(maxpowindx);
% % virtualgrid.leadfield = {[virtualgrid.leadfield{maxpowindx}]};
%
% cfg.method = 'pcc';
% cfg.pcc = cfg.dics;
% cfg = rmfield(cfg, 'dics');
% cfg.grid = virtualgrid;
% %cfg.rawtrial='yes';
% gammaChan = ft_sourceanalysis(cfg, freq_pst); % should the filter just be determined on (250ms or 750ms) gamma data?
% % here we have the single taper fourier coefficients in the mom-field
% mom = gammaChan.avg.mom{1};
% [u,s,v]=svd(real(mom*mom'));
% newmom=u(:,1)'*mom;
% newpow=abs(newmom).^2;
% newpow_trl=[];
% nTapers = size(mom,2)/length(ep_data.pst.trial);
% for iTrial = 1:length(ep_data.pst.trial)
%     newpow_trl(iTrial) = sum(newpow(iTrial*nTapers-(nTapers-1):iTrial*nTapers))/nTapers;
% end
% gammaPow = log(newpow_trl);
% gammaPow = (gammaPow-mean(gammaPow))/std(gammaPow);