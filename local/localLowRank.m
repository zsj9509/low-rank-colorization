function [resImg] = localLowRank( gImg, Obvs, mu, patPara)

patSize = patPara.patSize;
sliding = patPara.sliding;

% image to patch
imSize = size(gImg);
[pat, patIdx] = im2patch(gImg, patSize, sliding);
[patR, patG, patB] = colorIm2patch(Obvs, patSize, sliding);

% choose anchor point
anchor = sum((patR ~= -1), 1);
anchor = (anchor >= 1);
N = sum(anchor);

pat = pat(:,anchor);
patIdx = patIdx(anchor);

patR = patR(:,anchor);
patG = patG(:,anchor);
patB = patB(:,anchor);

% initial
resPatR = zeros(size(pat));
resPatG = zeros(size(pat));
resPatB = zeros(size(pat));

para.maxIter = 100;
para.tol = 1e-6;    
para.pnt = 0;

parfor n = 1:N  
    % optimization    
    gayPat = reshape(pat(:, n), patSize, patSize);
    gupObv = cat(2, reshape(patR(:,n), patSize, patSize), ...
        reshape(patG(:,n), patSize, patSize), ...
        reshape(patB(:,n), patSize, patSize));
    Omega = double(gupObv ~= -1);
    
    [clrPat, out] = optADMM( gayPat, gupObv, Omega, mu, para );
    
    clrPat = reshape(clrPat, patSize^2*3, 1);
    
    resPatR(:,n) = clrPat(1:patSize^2);
    resPatG(:,n) = clrPat(patSize^2+1:2*patSize^2);
    resPatB(:,n) = clrPat(2*patSize^2+1:end);

    fprintf('%d / %d, inner loop %d \n', n, N, length(out.obj));
end

resImg = patch2colorIm(resPatR, resPatG, resPatB, patIdx, patSize, imSize);

end