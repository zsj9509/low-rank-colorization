clear all; clc; close all;

mask = 0;
cImg = imread('images/castle.jpg');
cImg = double(cImg)/255;
[m, n, k] = size(cImg);

% observations
Omega = zeros(m, n);
idx = randperm(m*n);
idx = idx(1:floor(length(idx)*0.05));
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
figure;
imshow(gImg, []);

% obvs
Omega = repmat(Omega, 1, 3);
figure;
imshow(Omega, [])

clear temp1 temp2 idx;
close all;

% ------------------------------------------------------------------------

patPara.patSize = floor(sqrt(min(size(gImg)))/2);
patPara.kNN = min(floor(patPara.patSize^2/2), 50);
patPara.sliding = 2;
patPara.rho = 1;
patPara.pnt = 1;
patPara.epsilon = 2;

%-------------------------------------------------------------------------
propD = D;
propNum = 20;
nNnz = zeros(1, propNum);
propTol = 1e-3;
for i = 1:propNum
    [ Omega, propD ] = propColor(gImg, Omega, propD, propTol);
    nNnz(i) = nnz(Omega)/numel(Omega);
    if(i > 1 && nNnz(i) - nNnz(i - 1) < 0.005)
        break;
    end
    if(nnz(Omega)/numel(Omega) >= 0.1)
        break;
    end
    propTol = propTol/2;
end
clear propNum nNnz propTol;

for i = 1:10
    mu = 0.01*2^(i - 1);
%     [ rImg ] = localColorization( gImg, propD, mu, patPara);

    patPara.patSize = ceil(sqrt(min(size(gImg)))/2);
    patPara.sliding = 2;
    patPara.epsilon = min(patPara.patSize/3, 3);
    patPara.kNN = min(50, patPara.patSize^2);
    patPara.rho = 1;
    patPara.pnt = 1;

    rImg  = localLowRank(gImg, D, mu, patPara);

    PSNR(i) = psnr(rImg, cImg);
    SSIM(i) = ssim(rImg, cImg);
    recover{i} = rImg;

    save('castleMuLLR.mat');    
    
    fprintf('mu: %d \n', mu);
end

clear propD i;

