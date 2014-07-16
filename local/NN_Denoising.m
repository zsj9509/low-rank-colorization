function [ resImg, PSNR ] = Denoising( nosImg, clnImg, lambda )

patSize = 8;
sliding = 1;
imSize = size(nosImg);

[pat, patIdx, patPos] = im2patch(nosImg, patSize, sliding);

patNum = length(patIdx);

for repeat = 1:2
    if(repeat == 1)
        augPat = augPathce(pat, patPos, imSize, patSize/2);
    else
        resPat = 0.9*resPat + 0.1*pat;
        augPat = augPathce(resPat, patPos, imSize, patSize/2);
        lambda = 0.75*lambda;
    end
    kdTree = vl_kdtreebuild(augPat);
    resPat = zeros(size(pat));

    patWgt = zeros(1, patNum);
    for n = 1:patNum
        curPat = augPat(:, n);
        [gnpIdx, gnpDis] = vl_kdtreequery(kdTree, augPat, curPat, ...
            'NumNeighbors', patSize^2);

        gupPat = double(augPat(1:patSize^2, gnpIdx));
        meanPat = mean(gupPat, 2);
        meanPat = repmat(meanPat, 1, length(gnpIdx));

        gupPat = gupPat - meanPat;
        [newPat, patR] = patDenosing(gupPat, lambda);
        newPat = newPat + meanPat;
        newWgt = exp(-gnpDis');

        newPat = updateGupPat(newPat, resPat(:, gnpIdx), newWgt, patWgt(gnpIdx));
        resPat(:, gnpIdx) = newPat;

        patWgt(gnpIdx) = patWgt(gnpIdx) + newWgt;

        fprintf('%d out of %d, patch group rank %d  \n', n, patNum, patR);
    end
    resImg = patch2im(resPat, patIdx, patSize, imSize);
    PSNR = 20*log10(1/sqrt(mean((resImg(:) - clnImg(:)).^2)));
end

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
function [pat, r] = patDenosing(pat, lambda)

[U, S, V] = svd(pat, 'econ');

S = diag(S);
S = S - lambda;
S = S(S > 0);

r = length(S);

U = U(:, 1:r);
V = V(:, 1:r);
S = diag(S);

pat = U*S*V';

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
