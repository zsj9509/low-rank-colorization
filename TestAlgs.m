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
    idx = idx(1:floor(length(idx)*0.1));
    Omega(idx) = 1;
end
Omega = repmat(Omega, 1, 3);

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

clear temp1 temp2 idx;

% %% ------------------------------------------------------------------------

Data.Omega = Omega;
Data.D = reshape(double(D), m, n*k);
Data.B = gImg/3;

para.maxIter = 100000;
para.tol = 1e-8;

t = tic;
[ ~, obj_pro ] = optProximal( gImg, Data.D, Omega, 2, para );
t_pro = toc(t);

t = tic;
[ ~, obj_re  ] = optReweight( gImg, Data.D, Omega, 2, para );
t_re = toc(t);

t = tic;
[ ~, obj_alm ] = optADMM( gImg, Data.D, Omega, 2, para );
t_alm = toc(t);

close all;
figure;
p = loglog(obj_pro);
hold on;
set(p, 'color', 'red');
p = loglog(obj_re);
set(p, 'color', 'black');
p = loglog(obj_alm);
set(p, 'color', 'blue');
