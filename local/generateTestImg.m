function [gImg, D, Omega ] = generateTestImg( cImg, ratio )
% cImg: [m n 3] color image
% gImg: R + B + G
% D: labels [m n 3] position indicated by Omega [m n]

cImg = cImg + 1e-9;
[m, n, ~] = size(cImg);

% observations
Omega = zeros(m, n);
idx = randperm(m*n);
idx = idx(1:floor(length(idx)*ratio));
Omega(idx) = 1;

% labels
D = zeros(size(cImg));
temp1 = zeros(m, n);
temp2 = cImg(:,:,1);
temp1(idx) = temp2(idx);
D(:,:,1) = temp1;

temp1 = zeros(m, n);
temp2 = cImg(:,:,2);
temp1(idx) = temp2(idx);
D(:,:,2) = temp1;

temp1 = zeros(m, n);
temp2 = cImg(:,:,3);
temp1(idx) = temp2(idx);
D(:,:,3) = temp1;

% gray images
gImg = (cImg(:,:,1) + cImg(:,:,2) + cImg(:,:,3));

end

