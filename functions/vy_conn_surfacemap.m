function vy_conn_surfacemap(filenameSURF,FSfolder, connpath,spmpath, savename)

addpath(connpath);
addpath(spmpath);

h = get(0, 'Children');
if isempty(findobj(h,'tag','CONN functional connectivity toolbox'))
    conn
end

% filenameSURF = ['.\output\',savename,'.nii'];
filenameVOL = [];
% FSfolder = 'E:\My Matlab\Connectivity\Conn\conn17a\conn\utils\surf';
sphplots = [];
connplots = [];
facealpha = 1;
position = [-1 0  0];
defaultfilepath = ['.\output\',savename];
conn_mesh_display(filenameSURF,[],FSfolder)