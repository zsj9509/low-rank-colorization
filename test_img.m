clear all; clc; close all;
imName = {'castle', 'couple', 'koala', 'lake', 'landscape',...
    'mushroom', 'street', 'woman'};

for i = 1:1:8
for t = 1:1:1 % observations

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

% rImg = optADMM( Data.B*3, Data.D, Data.Omega, 0.01, para );
rImg = ColorizationLL(Data, 10, 10, 1e-2);
% rImg = ColorizationLR(Data, 10, 10);
imSize = size(cImg);
rImg = reshape(rImg, imSize(1), imSize(2), 3);
clear imSize;

Global_PSNR(i,t) = psnr(rImg, cImg);
recover{i,t}.Global = rImg;

%-------------------------------------------------------------------------
[m, n, k] = size(D);
propD = reshape(D, m, n*k);
propD = propD + randn(size(propD))*0.02;
propD = propD.*Omega;
propD = reshape(propD, m, n, k);

[ rImg ] = colorUseOpt(gImg/3, propD );
 
UseOpt_PSNR(i,t) = psnr(rImg, cImg);
recover{i,t}.UseOpt = rImg;

clear m n k propD;

%-------------------------------------------------------------------------
patPara.sliding = 2;
patPara.epsilon = 0.9;
patPara.rho = 1;
patPara.pnt = 1;

% propD = D;
% propNum = 10;
% nNnz = zeros(1, propNum);
% propTol = 1e-2/t;
% for i = 1:propNum
%     [ Omega, propD ] = propColor(gImg, Omega, propD, propTol);
%     nNnz(i) = nnz(Omega)/numel(Omega);
%     if(nnz(Omega)/numel(Omega) >= 0.10)
%         break;
%     end
%     propTol = propTol/2;
% end
% clear propNum nNnz propTol;

% -------------------------------------------------------------------------
patPara.patSize = 20;
nnzD = 0;
radius = 8;
while(nnzD < 0.03)
    radius = radius + 2;
    [propD] = LocalColorConsistency(Data.D, Data.B, Data.Omega, 1e-3, radius);
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

rImg  = localLowRank(gImg, propD, 0.04, patPara);
Local_PSNR(i,t) = psnr(rImg, cImg);
recover{i,t}.LLR = rImg;

% -------------------------------------------------------------------------
patPara.patSize = 16;
patPara.kNN = min(floor(patPara.patSize^2/4), 50);
[ rImg ] = localColorization( gImg, propD, (0.01*patPara.patSize), patPara);
clear propD m n Data;

GupLLR_PSNR(i,t) = psnr(rImg, cImg);
recover{i,t}.GupLLR = rImg;

% matName = strcat(imName, sprintf('-%.2f.mat', 0.01*t));
% save(matName);
save('16.mat');

end % repeat missing
end % image name


