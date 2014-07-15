clear all;

clnImg = imread('hotel.jpg');
[m, n, k] = size(clnImg);

patchSize = 8;
sliding = 1;

if(k > 1)
    clnImg = rgb2gray(clnImg);
end
clnImg = double(clnImg);
clnImg = clnImg/255;

[Pat,idx, pos] = im2patch(clnImg, patchSize, sliding);
c = patchSize/2;
pos(1,:) = c*pos(1,:)/n;
pos(2,:) = c*pos(2,:)/m;

Pat = single(Pat);
pos = single(pos);

AugPat = cat(1, Pat, pos);

kdTree = vl_kdtreebuild(AugPat);

N = ceil(rand()*length(idx));
[index, distance] = vl_kdtreequery(kdTree, AugPat, AugPat(:,N), 'NumNeighbors', patchSize^2) ;

close all;
figure;
PatGup = Pat(:,index);
I = patch2im(PatGup, idx(index), patchSize, [m, n]);

subplot(1, 2, 1);
imshow(uint8(255*I));

subplot(1, 2, 2);
s = svd(PatGup);
plot(s);