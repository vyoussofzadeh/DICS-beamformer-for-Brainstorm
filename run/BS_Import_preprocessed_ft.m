
%% Steps to import Preprocessed Fieldtrip import
% 1: Open bf
% 2: Manually import ft-ic-preprocessed data into BS (using: import EEG/MEG)
% 3: Export channel file (using individual bf datafile) into Matlab
% workspce. Export it as a new variable, e.g. Channels.
disp('These are completed using bs GUI');

%%
% 4: import ft-ic-preprocessed into matlab workspace. Run below script
clear cln_data
% ic_file = '/data/MEG/Vahab/test_data/processed/ft_process/older/raghavan_manoj2/dfn/ic_raghavan_manoj2';
ic_file = '/data/MEG/Clinical/ft_process/19/dougherty_danielle/DFN/ic_dougherty_danielle.mat';

load(ic_file);
disp('DONE');

%%
% 5: Editing channel file: run below script (while BS is open)
channels2 = channels;
channels2.Channel=channels2.Channel(1:306);

clear a b
a = cln_data.label;
for i=1:length(channels2.Channel)
    b{i,:} = channels2.Channel(i).Name;
end
[sharedvals,idx] = intersect(b,a,'stable');
channels2.Channel=channels2.Channel(idx);
channels2.Projector = [];
disp('DONE');

%%
% 6: Updatig channel file: import into bf (file/import from matlab and choose channels2)
% 7: run the averager (average files)
% 8: run/extract/extract time from -0.2 to 1.2sec
