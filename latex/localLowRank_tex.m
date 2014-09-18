clear all;

gImg = imread('images/lake.jpg');

if(length(size(gImg)) > 2)
    gImg = rgb2gray(gImg);
end

gImg = double(gImg)/255;

sliding = 1;
patSize = 10;
imSize = size(gImg);

[pat, patIdx, patPos] = im2patch(gImg, patSize, sliding);
pat = augPatch(pat, patPos, imSize, 2);
kdTree = vl_kdtreebuild(pat);

N = ceil(rand()*length(patIdx));
curPat = pat(:, N);
[gnpIdx, gnpDis] = vl_kdtreequery(kdTree, pat, curPat, ...
    'NumNeighbors', patSize^2/2);
gnpIdx = gnpIdx';
gupPat = pat(1:patSize^2, gnpIdx);

clear gnpIdx N;

close all;

sgbl = svd(gImg);
sgbl = sgbl(2:end);
sgbl = sgbl - min(sgbl);
slow = svd(double(gupPat));
slow = slow - min(slow);

curPat = curPat(1:patSize^2);
curPat = curPat + 0.1*randn(size(curPat));
curPat = reshape(curPat, patSize, patSize);
spat = svd(curPat);
spat = spat - min(spat);

maxS = max(sgbl(1), slow(1))*1.2;

%% ---------------------------------------------------------------
figure;
p = semilogy(slow + 1);
set(p, 'linewidth',2);
axis([1,min(size(gupPat)),0,maxS]);
grid on;
xlabel('The i-th singular value');
ylabel('Magnitude log(s + 1)');

figure;
p = semilogy(sgbl + 1);
set(p, 'linewidth',2);
axis([1,min(size(gImg)) - 1,0,maxS]);
grid on;
xlabel('The i-th singular value');
ylabel('Magnitude log(s + 1)');

figure;
p = semilogy(spat + 1);
set(p, 'linewidth',2);
axis([1,min(size(spat, 1)),0,maxS]);
grid on;
xlabel('The i-th singular value');
ylabel('Magnitude log(s + 1)');
