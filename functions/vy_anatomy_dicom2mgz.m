function vy_anatomy_dicom2mgz(subject)
%streams_anatomy_dicom2mgz takes the the subject info data structure (or subject string as 'sXX') 
%   
%   Picks up the dicom files, reslices the image and creates a .mgz file (spm coordsyst)

if ischar(subject)
  subject = streams_subjinfo(subject);
end

subject_code = subject.name;
subject_number = str2num(subject.name(2:end));
anatomy_savedir = fullfile('/project/3011044.02/preproc/anatomy');

% select the last dicom file in subject's mri directory
dicom_dir  = fullfile(subject.mridir);


dicom_subdir = dir(dicom_dir);
dicom_subdir = dicom_subdir(end).name;
dicom_list = dir(fullfile(dicom_dir, dicom_subdir)); % choose the second folder with anatomical dicoms
dicom_file = fullfile(dicom_dir, dicom_subdir, dicom_list(end).name);

% read in the dicom files
mri   = ft_read_mri(dicom_file);

% filename for saving
mgz_filename = fullfile(anatomy_savedir, [subject_code, '_mri' '.mgz']);

% save the images in the mgz format
cfg             = [];
cfg.filename    = mgz_filename;
cfg.filetype    = 'mgz';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri);

end

