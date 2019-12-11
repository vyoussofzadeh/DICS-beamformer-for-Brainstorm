disp('1: Parcellation, evc');
disp('2: Whole-brain, evc');
disp('3: Whole-brain, conn');

technique = input('Enter the method: ');

clear evc conn
k = 1;
for i=1:length(files_sel)
    disp(files_sel{i})
    d1 = dir(fullfile(DestDirectory,files_sel{i})); files = {d1.name}'; files_sel1 = files(3:end,1);
    for j=1:length(files_sel1)
        d2 = fullfile(DestDirectory,files_sel{i},files_sel1{j},'network');
        switch technique
            case 1
                files = dir(fullfile(d2,'network_evc_par*.mat'));
            case 2
                files = dir(fullfile(d2,'network_evc_C-*.mat'));
            case 3
                files = dir(fullfile(d2,'network_evc_conn*.mat'));
        end
        if isempty(files) == 0
            load(fullfile(files.folder,files.name));
            switch technique
                case 1
                    evc(k,:) = data_intpar.eigenvector_cent;
                case 2
                    evc(k,:) = network_diff_lcmv.eigenvector_cent;
                case 3
                    conn(k,:,:) = conn_ratio.plvspctrm ;
            end
            data_dis{k} = fullfile(files.folder,files.name);
            k=k+1;
            disp(k-1)
        end
    end
end

outputdir = fullfile(DestDirectory,'Group ave',mtd);
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end

evc_n = [];
switch technique
    case 1
        D = data_intpar;
        evcn = vy_normalize(evc); % zero-mean, divide by std, average
        D.eigenvector_cent = evcn';
        %         D.eigenvector_cent = (mean(evc,1)./std(evc))';
    case 2
        D = network_diff_lcmv;
        evcn = vy_normalize(evc);
        D.eigenvector_cent = evcn';
        %         D.eigenvector_cent = (mean(evc,1)./std(evc))';  
    case 3
        D = conn_ratio;
end

%%
switch technique
    case 1
        savepath = fullfile(outputdir,'eigenvec_group_par.mat');
        save(savepath, 'D','evc', '-v7.3');
    case 2
        savepath = fullfile(outputdir,'eigenvec_group_wb.mat');
        save(savepath, 'D','evc', '-v7.3');
    case 3
        savepath = fullfile(outputdir,'eigenvec_group_conn.mat');
end

%% Source vis - quick inspection
switch technique
    case {1,2}
        cfg              = [];
        cfg.method       = 'ortho';
        cfg.funparameter = 'eigenvector_cent';
        figure
        ft_sourceplot(cfg,D);
        hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
end

%%
% D.eigenvector_cent(D.eigenvector_cent < 0.7.*max(D.eigenvector_cent)) = NaN; % negative effects
% D.eigenvector_centdimord = 'chan';

%% Source vis - template
% D.eigenvector_cent = 100.*D.eigenvector_cent;
gtm = 'eigenvector_cent';
clear savepath
switch technique
    case 1
        savepath{1} = fullfile(outputdir,'n_par1_t_2_par');
        savepath{2} = fullfile(outputdir,'n_par2_t_3_par');
        savenii = fullfile(outputdir,'n_par.nii');
    case 2
        savepath{1} = fullfile(outputdir,'n_par1_t_2_wb');
        savepath{2} = fullfile(outputdir,'n_par2_t_3_wb');
        savenii = fullfile(outputdir,'n_wb.nii');
end
switch technique
    case {1,2}
        vy_mapvisualisation(D,gtm,0.5, savepath,0);
        vy_savenifti(D,gtm,savenii);
end

%%
switch technique
    case 1
        %- MNI coordinates
        coor = [];
        pp = D.brainordinate;
        for i=1:length(pp.tissuelabel)
            idx = pp.tissue == i;
            M = mean(pp.pos(idx,:),1);
            M = round(M*100)/100; % rounding into 2 decimals
            coor(i,:) = round(M);
        end
        [ROI, ROI_sel] = vy_ROI_report(D,.8, coor, gtm);
        disp(ROI_sel)
        savepath = fullfile(outputdir,'n_ROIs');
        hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
        
        textfile_rej = fullfile(outputdir,'ROI_sel');
        writetable(ROI_sel,textfile_rej,'Delimiter',' ');
        
    case 2
        gtm = 'eigenvector_cent';
        [~, D_par, coor] = vy_parcellate(D, atlas, gtm);
        D_par.eigenvector_centdimord = 'chan';
        
        savepath = fullfile(outputdir,'eigenvec_group_wb_par.mat');
        save(savepath, 'D_par', '-v7.3');
        
        clear savepath
        savepath{1} = fullfile(outputdir,'n_par1_t_2_wb_par');
        savepath{2} = fullfile(outputdir,'n_par2_t_3_wb_par');
        vy_mapvisualisation(D_par,gtm,0.7, savepath);
        
        savenii = fullfile(outputdir,'n_wb_par.nii');
        vy_savenifti(D_par, gtm, savenii);
        
        [ROI, ROI_sel] = vy_ROI_report(D_par,.8, coor, gtm);
        disp(ROI_sel)
        savepath = fullfile(outputdir,'n_ROIs_wb_par');
        hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
        
        textfile_rej = fullfile(outputdir,'ROI_sel_wb_par');
        writetable(ROI_sel,textfile_rej,'Delimiter',' ');
end

%% Individuals
switch technique
    case {1,2}
        switch technique
            case 1
                outputinddir = fullfile(outputdir,'indiv_n');
            case 2
                outputinddir = fullfile(outputdir,'indiv_n_wb');
        end
        if exist(outputinddir, 'file') == 0
            mkdir(outputinddir);   %create the directory
        end
        clear savenii
        for i=1:size(evc,1)
            
            Index = strfind(data_dis{i}, '\');
            d = data_dis{i};
            %     clear savepath
            %     savepath{1} = fullfile(outputinddir,['1_',d(Index(end)+1:end-4)]);
            %     savepath{2} = fullfile(outputinddir,['2_',d(Index(end)+1:end-4)]);
            D.eigenvector_cent = evc(i,:)';
            %     vy_mapvisualisation(D,gtm,0, savepath);
%             vy_mapvisualisation(D,gtm,0, []);
            
            savenii{i} = fullfile(outputinddir,[d(Index(end)+1:end-4),'.nii']);
            vy_savenifti(D,gtm,savenii{i});
            
            close all
        end
        save(['.\data_file\meg_',tsk],'savenii');
end

%% Conn
switch technique
    case 3
        %- conn normalize
        conn_n = conn;
        for i=1:size(conn,1)
            disp([num2str(i),'/',num2str(size(conn,1))])
            A = squeeze(conn_n(i,:,:));
            A = A./max(A(:));
            conn_n(i,:,:) = A;
        end
        
        %- Conn (avg)
        D.plvspctrm = squeeze(mean(conn_n,1));
        
        conn_par.conn_thre = 0.80;
        vy_connvis2(D,individual_headmodel,conn_par);
        view([156,47])
        
        savepath = fullfile(outputdir,'conn_avg');
        hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
        
        savepath = fullfile(outputdir,'conn_avg.mat');
        save(savepath, 'conn_ratio', 'individual_headmodel','-v7.3');
end