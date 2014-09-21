remain = 0.2;

img = double(img) / 256;
A = [img(:, :, 1), img(:, :, 2), img(:, :, 3)];

nrm_A = norm(A, 'fro');

idx = rand(481, 321);
idx(idx < remain) = 0;
idx = ~~idx;

Omega = [idx, idx, idx];
D = A;
D(Omega) = 0;

Red = D(:, 1: 321);
Green = D(:, 322: 642);
Blue = D(:, 643:963);

Data.Omega = idx;
Data.Red = Red; 
Data.Green = Green; 
Data.Blue = Blue; 
Data.BW = ( A(:, 1: 321) + A(:, 322: 642) + A(:, 643:963) ) / 3;

save Man.mat Data A nrm_A


% Show remaining entries
img = reshape([Data.Red, Data.Green, Data.Blue], [481, 321, 3]);imshow(img);

% compute
[R G B iter] = colorization(Data, 0.1, 0.1);

% restore image
img = reshape([R, G, B], [481, 321, 3]);imshow(img);

