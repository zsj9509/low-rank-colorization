clear all; clc; close all;

mask = 0;
cImg = imread('images/hotel-d.jpg');
cImg = double(cImg)/255;
[m, n, k] = size(cImg);

% observations
if(1 == mask)
    Omega = imread('images/hotel-d-mask2.bmp');
    Omega = double(Omega < 10);
    idx = find(Omega == 1);
else
    Omega = zeros(m, n);
    idx = randperm(m*n);
    idx = idx(1:floor(length(idx)*0.05));
    Omega(idx) = 1;
end
Omega = repmat(Omega, 1, 3);

figure;
imshow(Omega, [])

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

clear temp1 temp2 idx;
close all;

% %% ------------------------------------------------------------------------

% Data.Omega = Omega;
% Data.D = reshape(double(D), m, n*k);
% Data.B = gImg/3;
% 
% [rImg, ~, ~, iter] = ColorizationLR(Data, 10, 10, 1e-5, 500);
% 
% AAAI = norm(rImg(:) - cImg(:), 2)/norm(cImg(:), 2);
% figure;
% rImg = reshape(rImg, m, n, k);
% imshow(rImg, []);
% title('AAAI');

% para.maxIter = 10000;
% para.tol = 1e-10;
% para.pnt = 1;
% 
% [ rImg, obj_pro ] = optProximal( gImg, Data.D, Omega, 1, para );
% [ rImg, obj_re ] = optReweight( gImg, Data.D, Omega, 1, para );
% [ rImg, obj_alm ] = optADMM( gImg, Data.D, Omega, 1, para );
% 
% Global = norm(rImg(:) - cImg(:), 2)/norm(cImg(:), 2);
% figure;
% rImg = reshape(rImg, m, n, k);
% imshow(rImg, []);
% title('global');

[ rImg ] = localColorization( gImg, D, 0.3);

Local = norm(rImg(:) - cImg(:), 2)/norm(cImg(:), 2);
figure;
imshow(rImg, [])
title('local');

% [ rImg ] = colorUseOpt( gImg/3, D );
%  
% UseOpt = norm(rImg(:) - cImg(:), 2)/norm(cImg(:), 2);
% figure;
% imshow(rImg, [])
% title('useOpt');