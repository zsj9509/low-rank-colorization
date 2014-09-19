clear all; clc; close all;
imName = 'castle';

cImg = imread(strcat('images/', imName, '.jpg'));
cImg = double(cImg)/255;
[m, n, k] = size(cImg);

% observations
Omega = zeros(m, n);
idx = randperm(m*n);
idx = idx(1:floor(length(idx)*0.01));
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

gImg = (cImg(:,:,1) + cImg(:,:,2) + cImg(:,:,3));
Omega = repmat(Omega, 1, 3);

clear temp1 temp2 idx;
close all;

% ---------------------------------------------------------
Data.Omega = Omega;
Data.D = reshape(double(D), m, n*k);
Data.B = gImg/3;

% ---------------------------------------------------------
maxIter = 500;
tol = 1e-5;

% ---------------------------------------------------------
% l1
[G] = LocalColorConsistency(Data.D, Data.B, Data.Omega, 5e-2);
O = double(G > 0);
Data.Omega = O;
Data.D = G;

for eta = 1:15
    for lambda = 1:15
        eta_i = 0.01*2^eta/2;
        lambda_j = 0.01*2^lambda/2;
        
        rImg1 = ColorizationLR(Data, lambda_j, eta_i, tol);
        imSize = size(cImg);
        rImg1 = reshape(rImg1, imSize(1), imSize(2), 3);
        PSNR1(eta, lambda) = psnr(rImg1, cImg);
    end
    save;
end

clear eta lambda;

% l2
para.tol = tol;
para.maxIter = maxIter;
para.pnt = 1;
para.acc = 1;

for mu = 1:15
    for lambda = 1:15
        mu_i = 0.01*2^mu/2;
        para.lambda = 0.01*2^lambda/2;
        
        rImg2 = optADMM( Data.B*3, Data.D, Data.Omega, mu_i, para );
        imSize = size(cImg);
        rImg2 = reshape(rImg2, imSize(1), imSize(2), 3);
        PSNR2(mu, lambda) = psnr(rImg2, cImg);
    end
    save;
end

clear mu lambda;
clear eta_i mu_i lambda_j rImg2 rImg1 imSize;

save('l1_l2.mat');
