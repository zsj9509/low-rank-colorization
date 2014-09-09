clear all;

cImg = imread('images/hotel.jpg');
cImg = double(cImg)/255;
gImg = sum(cImg, 3)/3;

sliding = 2;
patSize = 10;
imSize = size(gImg);

[patR, patG, patB] = colorIm2patch(cImg, patSize, sliding);
[pat, patIdx, patPos] = im2patch(gImg, patSize, sliding);
pat = augPatch(pat/3, patPos, imSize, 2);
kdTree = vl_kdtreebuild(pat);

N = ceil(rand()*length(patIdx));
curPat = pat(:, N);
[gnpIdx, gnpDis] = vl_kdtreequery(kdTree, pat, curPat, ...
    'NumNeighbors', patSize^2);
gnpIdx = gnpIdx';

patR = patR(:,gnpIdx);
patG = patG(:,gnpIdx);
patB = patB(:,gnpIdx);
patIdx = patIdx(gnpIdx);

resImg = patch2colorIm(patR, patG, patB, patIdx, patSize, imSize);

close all;
figure;
imshow(resImg, []);
figure;
s = svd([patR, patG, patB]);
