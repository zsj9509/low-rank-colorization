clear all; clc; close all;
% imName = {'castle', 'couple', 'koala', 'lake', 'landscape',...
%     'mushroom', 'street', 'woman'};

imName = {'castle-s'};

for i = 1:1:1
for t = 9:1:9 % observations

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

% % rImg = optADMM( Data.B*3, Data.D, Data.Omega, 0.01, para );
% rImg = ColorizationLL(Data, 10, 10,0.9e-1);
% % rImg = ColorizationLR(Data, 10, 10);
% imSize = size(cImg);
% rImg = reshape(rImg, imSize(1), imSize(2), 3);
% clear imSize;
% 
% Global_PSNR(i,t) = psnr(rImg, cImg);
% % recover{i,t}.Global = rImg;

%-------------------------------------------------------------------------
% [m, n, k] = size(D);
% propD = reshape(D, m, n*k);
% propD = propD + randn(size(propD))*0.015;
% propD = propD.*Omega;
% propD = reshape(propD, m, n, k);
% 
% [ rImg ] = colorUseOpt(gImg/3, propD );
%  
% UseOpt_PSNR(i,t) = psnr(rImg, cImg);
% % recover{i,t}.UseOpt = rImg;
% 
% clear m n k propD;

%-------------------------------------------------------------------------
patPara.sliding = 2;
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
nnz(propD)/numel(propD)
propD(Omega == 0) = -1;
propD = reshape(propD, m, n, k);
clear m n k nnzD D Omega;

% rImg  = localLowRank(gImg, propD, 0.04, patPara);
% Local_PSNR(i,t) = psnr(rImg, cImg);
% recover{i,t}.LLR = rImg;

% -------------------------------------------------------------------------
% patPara.patSize = 19 - tempT;
patPara.patSize = 10;
if(i == 7)
    patPara.rand = 0;
    patPara.epsilon = 0.5;
else
    patPara.rand = 1;
    patPara.epsilon = 0.95;
end
patPara.kNN = min(floor(patPara.patSize^2/4), 50);
clear propD m n Data;

% [ rImg ] = localColorization( gImg, propD, (0.01*patPara.patSize), patPara);
% GupLLR_PSNR(i,t) = psnr(rImg, cImg);
% recover{i,t}.GupLLR = rImg;

[ rImg_l2 ] = localColorization( gImg, propD, (0.01*patPara.patSize), patPara);
[ rImg_l1 ] = localColorization_l1( gImg, propD, (0.01*patPara.patSize), patPara);

clear rImg;
% save('random.mat');

end % repeat missing
end % image name


