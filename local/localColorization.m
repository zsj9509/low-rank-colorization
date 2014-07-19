function [ resImg ] = localColorization( gImg, Obvs, mu )

patSize = 6;
sliding = 1;
imSize = size(gImg);

[pat, patIdx, patPos] = im2patch(gImg, patSize, sliding);
[patR, patG, patB] = colorIm2patch(Obvs, patSize, sliding);

patNum = length(patIdx);

pat = augPathce(pat, patPos, imSize, patSize/2);
kdTree = vl_kdtreebuild(pat);

resPatR = patR;
resPatG = patG;
resPatB = patB;

patWgt = zeros(1, size(pat, 2));
for n = 1:patNum
    % find neighbors
    curPat = pat(:, n);
    [gnpIdx, gnpDis] = vl_kdtreequery(kdTree, pat, curPat, ...
        'NumNeighbors', patSize^2);
    gnpIdx = gnpIdx';

    % opt
    gupGrayPat = double(pat(1:patSize^2, gnpIdx));
    gupObvsPat = [patR(:,gnpIdx), patG(:,gnpIdx), patB(:,gnpIdx)];

    Omega = double(gupObvsPat ~= -1);
    if(sum(Omega(:)) < ceil(patSize/2))
        continue;
    end

    if(exist('newPat','var'))
        [ newPat, gupRank, iter ] = optProximal( gupGrayPat, ...
            gupObvsPat, Omega, mu, newPat );
    else
        [ newPat, gupRank, iter ] = optProximal( gupGrayPat, ...
            gupObvsPat, Omega, mu );
    end
    iter = length(iter);

    newWgt = exp(-gnpDis');

    % update recover patches
    N = length(gnpIdx);
    resPatR(:, gnpIdx) = ...
        updateGupPat(newPat(:,1:N), resPatR(:, gnpIdx), newWgt, patWgt(gnpIdx));
    resPatG(:, gnpIdx) = ...
        updateGupPat(newPat(:,N+1:2*N), resPatG(:, gnpIdx), newWgt, patWgt(gnpIdx));
    resPatB(:, gnpIdx) = ...
        updateGupPat(newPat(:,2*N + 1:end), resPatB(:, gnpIdx), newWgt, patWgt(gnpIdx));

    patWgt(gnpIdx) = patWgt(gnpIdx) + newWgt;

    fprintf('%d out of %d, rank %d, obvs: %1.3f, iter: %d \n', ...
        n, patNum, gupRank, sum(Omega(:))/numel(Omega), iter);
end

resImgR = patch2im(resPatR, patIdx, patSize, imSize);
resImgG = patch2im(resPatG, patIdx, patSize, imSize);
resImgB = patch2im(resPatB, patIdx, patSize, imSize);

resImg = cat(3, resImgR, resImgG, resImgB);

end

%% --------------------------------------------------------------

%% --------------------------------------------------------------
function [rPat, gPat, bPat] = colorIm2patch(cImg, patSize, sliding)

R = cImg(:,:,1);
G = cImg(:,:,2);
B = cImg(:,:,3);

rPat  = im2patch(R, patSize, sliding);
gPat  = im2patch(G, patSize, sliding);
bPat  = im2patch(B, patSize, sliding);

end

%% --------------------------------------------------------------
function [AugPat] = augPathce(pat, pos, imSize, c)
    m = imSize(1);
    n = imSize(2);
    
    pos(1,:) = c*pos(1,:)/n;
    pos(2,:) = c*pos(2,:)/m;

    pat = single(pat);
    pos = single(pos);

    AugPat = cat(1, pat, pos);
end

%% --------------------------------------------------------------
function [pat] = updateGupPat(newPat, oldPat, newWgt, oldWgt)

dim = size(newPat, 1);

newWgt = repmat(newWgt, dim, 1);
oldWgt = repmat(oldWgt, dim, 1);

newPat = newPat.*newWgt;
oldPat = oldWgt.*oldPat;

pat = newPat + oldPat;
pat = pat ./ (newWgt + oldWgt);

end

