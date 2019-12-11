function LI = vy_laterality(data)

% idx = [4:2:20,24:2:26,80:2:90];
% idx = [2:2:90];
% idx = [12:2:20,80:2:88];

idx = [12:2:20,80:2:88];
% idx(6) = []; idx(5) = [];
% idx = [idx-1, idx];

clear rightFT
for i=1:length(idx)
    rightFT{i,:} = data.label{idx(i)};
end
disp('left-side ROIs:')
disp(rightFT)
m_right  = mean(data.value(:,idx),2);

% Left FT lobe
idx = idx-1;
% idx = [1:2:90];

clear leftFT
for i=1:length(idx)
    leftFT{i,:} = data.label{idx(i)};
end
disp('left-side ROIs:')
disp(leftFT)

m_left  = mean(data.value(:,idx),2);

LI = (m_left - m_right)./ (m_left + m_right);
disp('Laterality per samples:')
disp(LI)




