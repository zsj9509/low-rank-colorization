bw = imread('BW.png');
label = imread('manual.png');

[m n p] = size(bw);

average = double(bw) / 256;
average = (average(:, :, 1) + average(:, :, 2) + average(:, :, 3)) / 3;

bw = [bw(:, :, 1), bw(:, :, 2), bw(:, :, 3)];
label = [label(:, :, 1), label(:, :, 2), label(:, :, 3)];
Omega = (abs(bw - label) > 0.00001);
bw = double(bw) / 256;
label = double(label) / 256;

% which pixels are labeled
idx = ~~(Omega(:, 1:n) + Omega(:, n+1:2*n) + Omega(:, 2*n+1:3*n));

% fit the label to the monochrome image
label_new = ( label(:, 1:n) + label(:, n+1:2*n) + label(:, 2*n+1:3*n) ) / 3;
normalizer = (average + 1e-8) ./ (label_new + 1e-8);
D = label .* [normalizer, normalizer, normalizer];

Data.D = D;
Data.Omega = [idx, idx, idx];
Data.B = average;

save manual.mat Data