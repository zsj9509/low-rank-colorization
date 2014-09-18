clear all; clc; close all;
imName = 'mushroom';

cImg = imread(strcat('images/', imName, '.jpg'));
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
% figure;
% imshow(gImg, []);

% obvs
Omega = repmat(Omega, 1, 3);
% figure;
% imshow(Omega, [])

clear temp1 temp2 idx;
close all;

Data.Omega = Omega;
Data.D = reshape(double(D), m, n*k);
Data.B = gImg/3;

% ------------------------------------------------------------------------
[L, S,~,~,O] = ColorizationLL(Data, 10, 10, 2e-2);
imSize = size(cImg);
rImg = L;
rImg = reshape(rImg, imSize(1), imSize(2), 3);



