function [ resImg ] = localColorization( gImg, Obvs, mu )

patSize = 7;
sliding = 1;
imSize = size(gImg);

[pat, patIdx, patPos] = im2patch(gImg, patSize, sliding);
[patR, patG, patB] = colorIm2patch(Obvs, patSize, sliding);

patNum = length(patIdx);

pat = augPatch(pat/3, patPos, imSize, 2);
kdTree = vl_kdtreebuild(pat);

resPatR = patR;
resPatG = patG;
resPatB = patB;

patWgt = zeros(1, size(pat, 2));
patJmp = 100;
for n = 1:patJmp:patNum
    patJmp = min(patJmp, patNum - n);
    % index and weight
    gnpIdx = zeros(patSize^2, patJmp);
    gnpDis = zeros(patSize^2, patJmp);
    curPat = pat(:, n:n + patJmp);
    parfor m = 1:patJmp
        curPat_m = curPat(:, m);
        [gnpIdx_m, gnpDis_m] = vl_kdtreequery(kdTree, pat, curPat_m, ...
            'NumNeighbors', patSize^2);
        gnpIdx_m = gnpIdx_m';
        
        gnpIdx(:, m) = gnpIdx_m;
        gnpDis(:, m) = exp(-2*gnpDis_m');
    end
    clear curPat curPat_m gnpIdx_m gnpDis_m;
    
    % observations
    gupObv = cell(1, patJmp);
    gayPat = cell(1, patJmp);
    for m = 1:patJmp
        gnpIdx_m = gnpIdx(:, m);
        
        gupObv_m = [patR(:,gnpIdx_m), patG(:,gnpIdx_m), patB(:,gnpIdx_m)];
        gayPat_m = 3*double(pat(1:patSize^2, gnpIdx_m));
        
        gupObv{m} = gupObv_m;        
        gayPat{m} = gayPat_m;
    end
    clear gupObv_m gnpIdx_m gayPat_m;
    
    % optimization
    para.maxIter = 100;
    para.tol = 1e-5;    
    para.pnt = 0;
    
    newPat = cell(1, patJmp);
    iter   = zeros(1, patJmp); 
    parfor m = 1:patJmp
        Omega_m = double(gupObv{m} ~= -1);
        if(sum(Omega_m(:)) < ceil(patSize/2))
            continue;
        end
        Omega_m = sparse(Omega_m);
        gayPat_m = gayPat{m};
        gupObv_m = gupObv{m};

        [newPat{m}, iter_m] = optADMM( gayPat_m, gupObv_m, ...
            Omega_m, mu, para );
        iter(m) = length(iter_m);
    end
    iter = mean(iter);
    clear Omega_m gayPat_m gupObv_m;
    
    % weighted update results
    for m = 1:patJmp
        gnpIdx_m = gnpIdx(:, m);
        gnpDis_m = gnpDis(:, m);
        newPat_m = newPat{m};
        
        % update recover patches
        N = length(gnpIdx_m);
        resPatR(:, gnpIdx_m) = ...
            updateGupPat(newPat_m(:,1:N), resPatR(:, gnpIdx_m), gnpDis_m', patWgt(gnpIdx_m));
        resPatG(:, gnpIdx_m) = ...
            updateGupPat(newPat_m(:,N+1:2*N), resPatG(:, gnpIdx_m), gnpDis_m', patWgt(gnpIdx_m));
        resPatB(:, gnpIdx_m) = ...
            updateGupPat(newPat_m(:,2*N + 1:end), resPatB(:, gnpIdx_m), gnpDis_m', patWgt(gnpIdx_m));

        patWgt(gnpIdx_m) = patWgt(gnpIdx_m) + gnpDis_m';
    end

    fprintf('%d out of %d, avg loops in opt %2.2f \n', n, patNum, iter);
end

resImgR = patch2im(resPatR, patIdx, patSize, imSize);
resImgG = patch2im(resPatG, patIdx, patSize, imSize);
resImgB = patch2im(resPatB, patIdx, patSize, imSize);

resImg = cat(3, resImgR, resImgG, resImgB);

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

