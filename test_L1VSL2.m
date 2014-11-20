clear; clc; close all;
imName = {'castle', 'couple', 'koala', 'lake', 'landscape',...
    'mushroom', 'street', 'woman'};

for i = 1:length(imName)
for t = 10:10 % observations

imName_i = imName{i};
cImg = imread(strcat('images/', imName_i, '.jpg'));
cImg = double(cImg)/255;
[gImg, D, Omega ] = generateTestImg( cImg, 0.01*t );
Omega = repmat(Omega, 1, 3);

% ------------------------------------------------------------------------
Data.Omega = Omega;
[m, n, k] = size(D);
Data.D = reshape(double(D), m, n*k);
Data.B = gImg/3;
clear m n k;

%-------------------------------------------------------------------------
patPara.sliding = 3;
patPara.rho = 1;
patPara.pnt = 1;

tempT = floor(sqrt(t));
patPara.patSize = 20;
nnzD = 0;
radius = 9 - tempT;
while(nnzD < 0.03)
    radius = radius + 2;
    [propD] = LocalColorConsistency(Data.D, Data.B, Data.Omega, ...
        1e-3/tempT, radius);
    nnzD = nnz(propD)/numel(propD);
end
[m, n, k] = size(D);
D = reshape(D, m, n*k);
Omega = Omega + double(propD > 0);
propD = propD + D;
fprintf('nz: %d \n', nnz(propD)/numel(propD));
propD(Omega == 0) = -1;
propD = reshape(propD, m, n, k);
clear m n k nnzD D Omega;

% -------------------------------------------------------------------------
patPara.patSize = 19 - tempT;
patPara.rand = 0;
patPara.epsilon = 0.25;
% if(i == 7)
%     patPara.rand = 0;
%     patPara.epsilon = 0.25;
% else
%     patPara.rand = 1;
%     patPara.epsilon = 0.5;
% end
patPara.kNN = min(floor(patPara.patSize^2/4), 50);e3
clear m n Data tempT;

time = tic;
[ rImg ] = localColorization( gImg, propD, (0.01*patPara.patSize), patPara);
TIMEL2(i,t) = toc(time);
PSNRL2(i,t) = psnr(rImg, cImg);

time = tic;
[ rImg ] = localColorization_l1( gImg, propD, (0.01*patPara.patSize), patPara);
TIMEL1(i,t) = toc(time);
PSNRL1(i,t) = psnr(rImg, cImg);

clear rImg Omega cImg gImg D Data time;
save('time.mat');

end % repeat missing
end % image name


