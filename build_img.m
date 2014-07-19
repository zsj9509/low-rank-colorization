clear all; clc; close all;

cImg = imread('images/hotel-ss.jpg');
cImg = double(cImg)/255;
[m, n, k] = size(cImg);

% observations
idx = randperm(m*n);
idx = idx(1:floor(length(idx)*0.1));

Omega = zeros(m, n);
Omega(idx) = 1;
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

figure;
imshow(uint8(D))

% gray images
gImg = (cImg(:,:,1) + cImg(:,:,2) + cImg(:,:,3));
figure;
imshow(gImg, []);

Data.Omega = Omega;
Data.D = reshape(double(D), m, n*k);
Data.B = gImg/3;

clear temp1 temp2 idx;
close all;

[rImg, ~, ~, iter] = ColorizationLR(Data, 10, 10, 1e-5, 500);

AAAI = norm(rImg(:) - cImg(:), 2)/norm(cImg(:), 2);
figure;
rImg = reshape(rImg, m, n, k);
imshow(rImg, []);
title('AAAI');

[ L ] = optProximal( gImg, Data.D, Omega, 1 );

Global = norm(rImg(:) - cImg(:), 2)/norm(cImg(:), 2);
figure;
rImg = reshape(L, m, n, k);
imshow(rImg, []);
title('global');

[ rImg ] = localColorization( gImg, D, 0.3);

Local = norm(rImg(:) - cImg(:), 2)/norm(cImg(:), 2);
figure;
imshow(rImg, [])
title('local');

[ rImg ] = colorUseOpt( gImg/3, D );

UseOpt = norm(rImg(:) - cImg(:), 2)/norm(cImg(:), 2);
figure;
imshow(rImg, [])
title('useOpt');