% data = s_data.bsl;
% data = s_data.pst;
function stat = vy_plv(cfg, data)

data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq' 'source'});
inparam = 'crsspctrm';
outparam = 'plvspctrm';
normrpt = 1;

dtype = ft_datatype(data);

cfg.refindx     = ft_getopt(cfg, 'refindx', 'all');
cfg.refindx     = 1:size(data.pos,1);
[data, powindx, hasrpt] = univariate2bivariate(data, 'mom', 'crsspctrm', dtype, 'cmb', cfg.refindx, 'keeprpt', 0);

cfg.pchanindx = [];
cfg.allchanindx = [];

data.(inparam) = reshape(data.(inparam), [1 size(data.(inparam))]);

data.dimord = ['rpt_', data.dimord];
data.dimord = 'rpt_pos_pos_freq';
cfg.complex = 'abs';
% phase locking value

hasjack = 0;
normpow = 1;
cfg.feedback = 'none';
optarg = {'complex', cfg.complex, 'dimord', data.dimord, 'feedback', cfg.feedback, 'pownorm', normpow, 'hasjack', hasjack};
if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
[datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});


stat = keepfields(data, {'pos', 'dim', 'transform', 'inside', 'outside'});
stat.(outparam) = datout;
if ~isempty(varout)
    stat.([outparam, 'sem']) = (varout/nrpt).^0.5;
end

% just copy them over, alhtough we don't know for sure whether they are needed in the output
if isfield(data, 'freq'), stat.freq = data.freq; end
if isfield(data, 'time'), stat.time = data.time; end

if isfield(data, 'grad'), stat.grad = data.grad; end
if isfield(data, 'elec'), stat.elec = data.elec; end
if exist('dof', 'var'), stat.dof = dof; end

