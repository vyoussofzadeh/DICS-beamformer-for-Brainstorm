function data_fix = vy_ica_cleaning(f_data, lay, savepath,saveflag)

n = 20; % ICs

satis = 0;
disp('ica cleaning ...');
if exist(savepath, 'file') == 2
    load(savepath)
else
    cfg = [];
    cfg.metric = 'zvalue';  % use by default zvalue method
    cfg.latency = [-400,900];
    cfg.layout   = lay;   % this allows for plotting individual trials
    r_data   = ft_rejectvisual(cfg, f_data);
    
    comp = vy_ica(r_data,lay, n);
    title(savepath)
    
    %% ECG identifer
    %     ecg1 = ecg_data.trial{1, 1};
    %     clear r
    %     for i = 1:n
    %         ic1(i,:) = comp.trial{1, 4}(i,:);
    %         r(i) = corr2(ic1(i,:),abs(ecg1));
    %         [a,b]=sort(r)
    %     end
    %     figure, plot(ic1(14,:))
    %     figure, plot(ic1(4,:))
    %     figure, plot(ecg1)
    
    %%
    
    %rej component
    cfg = [];
    cfg.updatesens = 'no';
    bic = input('Select bad ICs for slected data:');
    cfg.component = comp.label(bic);
    data_fix = ft_rejectcomponent(cfg, comp, r_data);
    while (satis == 0)&&(isempty(bic) == 0)
        close all
        comp = vy_ica(data_fix,lay,n);
        satis = input('Statisfied with ICA, yes = 1, not yet = 0:');
        if satis == 1, close('all'), break, end
        bic = input('Select bad ICs for slected data:');
        cfg.component = comp.label(bic);
        data_fix = ft_rejectcomponent(cfg, comp, data_fix);
    end
    close all
    if saveflag ==1
        save(savepath, 'data_fix', '-v7.3');
    end
end