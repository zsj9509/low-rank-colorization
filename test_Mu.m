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
Omega = repmat(Omega, 1, 3);

clear temp1 temp2 idx;

%-------------------------------------------------------------------------
Data.Omega = Omega;
Data.D = reshape(double(D), m, n*k);
Data.B = gImg/3;
[propD] = LocalColorConsistency(Data.D, Data.B, Data.Omega, 1e-3, 10);
[m, n] = size(gImg);
propD = reshape(propD, m, n, 3);
nnz(propD)/numel(propD)
propD(propD == 0) = -1;
clear Data m n;

for i = 1:7
    mu = (5e-4)*10^(i - 1);
    
    patPara.patSize = 12;
    patPara.sliding = 3;
    patPara.epsilon = 0.5;
    patPara.kNN = 30;
    patPara.rho = 1;
    patPara.pnt = 1;

    % rImg  = localLowRank(gImg, propD, mu, patPara);
    [ rImg ] = localColorization( gImg, propD, mu, patPara);

    PSNR(i) = psnr(rImg, cImg);
    recover{i} = rImg;

    save('castleMuLLR.mat');    
    
    fprintf('mu: %d \n', mu);
end

clear propD i;

