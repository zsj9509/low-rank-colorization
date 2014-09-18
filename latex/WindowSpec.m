clear all;

img = imread('images/window.png');
img = double(img)/255;

[m, n, k] = size(img);
img = reshape(img, m, n*k);

[U, s, V] = svd(img, 'econ');
s = diag(s) + 0.00001;
% s = s - min(s);

minSize = length(s);
s = s';
s = s(4:end);

figure;
p = semilogy(s);
set(p, 'linewidth',2);
axis([1,minSize,0, max(s)]);
grid on;
xlabel('The i-th singular value');
ylabel('Magnitude log(s + 1)');