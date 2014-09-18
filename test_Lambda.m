clear all; clc; close all;

mask = 0;
cImg = imread('images/castle.jpg');
cImg = double(cImg)/255;
[m, n, k] = size(cImg);

% observations
Omega = zeros(m, n);
idx = randperm(m*n);
idx = idx(1:floor(length(idx)*0.10));
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

%-------------------------------------------------------------------------
propD = D;
propNum = 20;
nNnz = zeros(1, propNum);
propTol = 1e-4;
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

for i = 11:15
    patPara.lambda = 0.01*2^(i - 1);
    patPara.patSize = 12;
    patPara.sliding = 2;
    patPara.epsilon = 0.5;
    patPara.kNN = 30;
    patPara.rho = 1;
    patPara.pnt = 1;
    
    [ rImg ] = localColorization( gImg, propD, 0.16, patPara);

    %rImg  = localLowRank(gImg, D, mu, patPara);

    PSNR(i) = psnr(rImg, cImg);
    recover{i} = rImg;

    save('castleLambda.mat');    
    
    fprintf('lambda: %d \n', patPara.lambda);
end

clear propD i;

