function xn = vy_normalize(x)

if length(size(x))>3
    SEM = squeeze(std(x(),1))./sqrt(size(x,1));
    Mean = squeeze(mean(x(),1));
    xn = Mean./SEM;
    
else
%     SD = squeeze(std(x(),1));
    SEM = std(x)./sqrt(length(x));
    Mean = squeeze(mean(x(),1));
    xn = Mean./SEM;
end
% for i=1:size(x,1)
%     x1(i,:) = (x(i,:) - mean(x(i,:)))./std(x(i,:));
% end
% xn = (mean(x1,1))';

end