
clear; clc, close('all'); warning off

brainstorm
addpath(genpath('./functions'));

%%
datadir = '/data/MEG/Clinical/MEG/**/brainstorm_db/anat/**/subjectimage*.mat';
d = rdir(datadir);

clear datafolder datafile
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
end
datafolder_all = datafolder';
disp(datafolder_all)

datadir = '/data/MEG/Clinical/MEG/**/brainstorm_db/anat/@default_subject/subjectimage*.mat';
d = rdir(datadir);

clear datafolder datafile_ex
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile_ex{i} = d(i).name;
end
datafolder_ex = datafolder';
disp(datafolder_ex)

%
Acommon = intersect(datafolder_all,datafolder_ex);
datafolder_sel = setxor(datafolder_all,Acommon);
disp(datafolder_sel)

%%
for i = 1:size(datafolder_sel,1)
    datafile = datafolder_sel{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    disp(datafile)
    d = rdir(fullfile(datafile,'subjectimage*.mat'));
    if ~isempty(d)
        sMRI = d.name;
        cd(datafile)
        cd ..
        OutputFile = fullfile(pwd,'T1.nii');
        out_mri_nii( sMRI, OutputFile);
    end
end