clear all; clc; close all;

imName_i = 'castle';
cImg = imread(strcat('images/', imName_i, '.jpg'));
cImg = double(cImg)/255;
[gImg, D, Omega ] = generateTestImg( cImg, 0.05 );
Omega = repmat(Omega, 1, 3);

% ------------------------------------------------------------------------
Data.Omega = Omega;
[m, n, k] = size(D);
Data.D = reshape(double(D), m, n*k);
Data.B = gImg/3;
clear m n k;

% t = tic;
% ColorizationLL(Data, 10, 10,0.9e-1);
% Global_T = toc(t)*15;
% clear imSize;

%-------------------------------------------------------------------------
patPara.sliding = 2;
patPara.rho = 1;
patPara.pnt = 1;

% -------------------------------------------------------------------------
tempT = floor(sqrt(5));
patPara.patSize = 14;
nnzD = 0;
radius = 9 - tempT;
while(nnzD < 0.03)
    radius = radius + 2;
    [propD] = LocalColorConsistency(Data.D, Data.B, Data.Omega, ...
        1e-3/tempT, radius - tempT);
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

% t = tic;
% localLowRank(gImg, propD, 0.04, patPara);
% Local_T = toc(t);

% -------------------------------------------------------------------------
patPara.rand = 1;
patPara.epsilon = 0.5;

patPara.kNN = min(floor(patPara.patSize^2/4), 50);
t = tic;
localColorization( gImg, propD, (0.01*patPara.patSize), patPara);
GupLLR = toc(t);

clear t tempT radius propD ans;


