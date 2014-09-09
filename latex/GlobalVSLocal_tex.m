clear all;

img = imread('images/barbara.png');

img = double(img)/255;

[m, n, k] = size(img);
r = floor(min(m,n)*0.05);

col = reshape(img, m, n*k);

[U, S, V] = svd(col);

U = U(:, 1:r);
V = V(:, 1:r);
S = S(1:r, 1:r);
col = U*S*V';

col = reshape(col, m, n, k);

close all;
figure;
imshow(img, []);
title('original');
figure;
imshow(col, []);
title('low rank');