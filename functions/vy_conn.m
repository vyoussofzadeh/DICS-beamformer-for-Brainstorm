function [source_conn, network, gtm] = vy_conn(source,conn_par,net_par)

cfg         = [];
cfg.method = conn_par.method;
idx = conn_par.idx;
cfg.complex = conn_par.complex;
source_conn = ft_connectivityanalysis(cfg, source);
% source_conn1 = vy_plv(cfg, source);
source_conn.dimord    = 'pos_pos';

%-
gtm = net_par.gtm;

cfg = [];
cfg.method    = gtm;
cfg.parameter = idx;
cfg.threshold = net_par.threshold;
network = ft_networkanalysis(cfg,source_conn);
network.pos     = source_conn.pos;
% network.dim     = source_conn.dim;
network.inside  = source_conn.inside;