clear all; clc; close all;
imName = 'castle-s';

cImg = imread(strcat('images/', imName, '.png'));
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
maxIter = 10000;
tol = 1e-8;
% -l2 unacc
para.tol = tol;
para.maxIter = maxIter;
mu = 10;
para.lambda = 10;
para.pnt = 1;

para.acc = 0;
t0 = tic;
[ ~, output0 ] = optADMM( Data.B*3, Data.D, Data.Omega, mu, para );
t0 = toc(t0);

% -l2 acc
para.acc = 1;
t1 = tic;
[ ~, output1 ] = optADMM( Data.B*3, Data.D, Data.Omega, mu, para );
t1 = toc(t1);


