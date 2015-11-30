clear; clc; close all;
imName = 'castle';

for t = 5:1:5 % observations

cImg = imread(strcat('images/', imName, '.jpg'));
cImg = double(cImg)/255;
[m, n, k] = size(cImg);

% observations
Omega = zeros(m, n);
idx = randperm(m*n);
idx = idx(1:floor(length(idx)*0.01*t));
Omega(idx) = 1;

% labels
D = zeros(size(cImg));
temp1 = -ones(m, n);
temp2 = cImg(:,:,1);
temp1(idx) = temp2(idx);
D(:,:,1) = temp1;

temp1 = -ones(m, n);
temp2 = cImg(:,:,2);
temp1(idx) = temp2(idx);
D(:,:,2) = temp1;

temp1 = -ones(m, n);
temp2 = cImg(:,:,3);
temp1(idx) = temp2(idx);
D(:,:,3) = temp1;

% gray images
gImg = (cImg(:,:,1) + cImg(:,:,2) + cImg(:,:,3));
% figure;
% imshow(gImg, []);

% obvs
Omega = repmat(Omega, 1, 3);
% figure;
% imshow(Omega, [])

clear temp1 temp2 idx;
close all;

% % ------------------------------------------------------------------------
% Data.Omega = Omega;
% Data.D = reshape(double(D), m, n*k);
% Data.B = gImg/3;
% 
% rImg = ColorizationLL(Data, 5, 10, 5e-2);
% imSize = size(cImg);
% rImg = reshape(rImg, imSize(1), imSize(2), 3);
% 
% Global_PSNR = psnr(rImg, cImg);
% Global_SSIM = ssim(rImg, cImg);
% recover.Global = rImg;

% %-------------------------------------------------------------------------
% [ rImg ] = colorUseOpt(gImg/3, D );
%  
% UseOpt_PSNR = psnr(rImg, cImg);
% UseOpt_SSIM = ssim(rImg, cImg);
% recover.UseOpt = rImg;

%-------------------------------------------------------------------------
patPara.patSize = ceil(sqrt(min(size(gImg)))/1.5);
patPara.sliding = 2;
patPara.epsilon = 0.5;
patPara.kNN = min(50, patPara.patSize^2);
patPara.rho = 1;
patPara.pnt = 1;
patPara.cImg = cImg;

% rImg  = localLowRank(gImg, D, 0.04, patPara);
% 
% Local_PSNR = psnr(rImg, cImg);
% Local_SSIM = ssim(rImg, cImg);
% recover.LLR = rImg;

% -------------------------------------------------------------------------
propD = D;
propNum = 10;
nNnz = zeros(1, propNum);
propTol = 1e-2/t;
for i = 1:propNum
    [ Omega, propD ] = propColor(gImg, Omega, propD, propTol);
    nNnz(i) = nnz(Omega)/numel(Omega);
    if(nnz(Omega)/numel(Omega) >= 0.10)
        break;
    end
    propTol = propTol/2;
end
clear propNum nNnz propTol;

[ rImg ] = localColorization( gImg, propD, 0.16, patPara);

clear propD i m n Data;

GupLLR_PSNR = psnr(rImg, cImg);
GupLLR_SSIM = ssim(rImg, cImg);
recover.GupLLR = rImg;

matName = strcat(imName, sprintf('-%.2f.mat', 0.01*t));
save(matName);

end % repeat missing


