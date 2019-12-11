function trl = streams_definetrial(dataset, name)

[status, filename] = streams_existfile([name,'_trl.mat']);
if status
  load(filename);
else
  fprintf('computing trial definition for subject %s\n', name);
  
  % now we make a distinction between s10 or less, and s11 or more
  
  subjnr = str2double(name(2:end));
  if subjnr<=10,
    % for the subjects recorded in 2013
    if ~iscell(dataset)
      % convert to cell, to accommodate for the fact that some sessions may
      % consist of more than one dataset
      dataset = {dataset};
    end
    
    for kk = 1:numel(dataset)
      event = ft_read_event(dataset{kk});
      type  = {event.type}';
      sel   = strcmp(type, 'UPPT001');
      event = event(sel);
      
      val   = [event.value]';
      smp   = [event.sample]';
      
      trl = zeros(0,4);
      for k = 1:numel(val)-2
        tmp1 = val(k);
        tmp2 = val(k+1);
        tmp3 = val(k+2);
        
        % for a triplet of 1's, the first two are assumed to be the begin and
        % end triggers (provided this has been the first stimulus in the
        % experiment it's true)
        if all([tmp1 tmp2 tmp3]==1)
          begsmp = smp(k);
          endsmp = smp(k+1);
          offset = 0;
          tmptrl = [begsmp endsmp offset val(k)];
        elseif tmp1==1 && tmp2==tmp3 && tmp2~=1
          begsmp = smp(k+1);
          endsmp = smp(k+2);
          offset = 0;
          tmptrl = [begsmp endsmp offset val(k+1)];
        else
          % skip
          continue;
        end
        trl = [trl; tmptrl];
      end
      alltrl{kk} = trl;
    end
    trl = alltrl;
    
    if numel(dataset)==1
      trl = trl{1};
    end
  else
    
    % for the subjects recorded in 2016
    if ~iscell(dataset)
      % convert to cell, to accommodate for the fact that some sessions may
      % consist of more than one dataset
      dataset = {dataset};
    end
    
    for kk = 1:numel(dataset)
      event = ft_read_event(dataset{kk});
      type  = {event.type}';
      sel   = strcmp(type, 'UPPT001');
      event = event(sel);
      
      val   = [event.value]';
      smp   = [event.sample]';
      
      selbeg = find(ismember(val, [10 20 30 40 50]));
      selend = find(ismember(val, [10 20 30 40 50]+7));
      
      trl = [smp(selbeg) smp(selend) zeros(numel(selbeg),1) val(selbeg)];
      
      alltrl{kk} = trl;
    end
    trl = alltrl;
    
    if numel(dataset)==1
      trl = trl{1};
    end
  end
  
  filename = fullfile('/project/3011044.02/preproc/meg', [name,'_trl.mat']);
  save(filename, 'trl');
end
