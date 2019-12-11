function [int_data,badchannels] = vy_interpolate_meg(r_data, raw_data, neighbours, plotting)

badchannels = setdiff(raw_data.label, r_data.label);

% load('bti248_neighb');

if plotting == 1
    % plotting neighbours for inspection
    cfg            = [];
    cfg.grad = raw_data.grad;
    cfg.neighbours = neighbours;
    % cfg.senstype   = 'MEG';
    ft_neighbourplot(cfg, r_data);
end

if ~isempty(badchannels)
    % Interpolate channels
    cfg = [];
    cfg.grad                      = raw_data.grad;
    cfg.method                    = 'spline';
    cfg.neighbours                = neighbours;
    cfg.badchannel                = badchannels;
    % cfg.senstype                  = 'MEG';
    int_data = ft_channelrepair(cfg, r_data);
    
    if plotting == 1
        %- data inspection
        cfg = [];
        cfg.viewmode = 'vertical';
        cfg.continuous = 'no';
        cfg.channel = badchannels;
        ft_databrowser(cfg,int_data);
    end
    
else
    disp('no electrode was interpolated');
    int_data = r_data;
end

end
