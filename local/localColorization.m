function [ resImg ] = localColorization( gImg, Obvs, mu )

patSize = 14;
kNN = min(100, patSize^2);
sliding = 2;
imSize = size(gImg);

[pat, patIdx, patPos] = im2patch(gImg, patSize, sliding);
[patR, patG, patB] = colorIm2patch(Obvs, patSize, sliding);

patNum = length(patIdx);

pat = augPatch(pat/3, patPos, imSize, 2);
kdTree = vl_kdtreebuild(pat);

resPatR = zeros(size(patR));
resPatG = zeros(size(patG));
resPatB = zeros(size(patB));

patWgt = zeros(1, size(pat, 2));
patJmp = min(floor(patNum*0.05), 400);
for n = 1:patJmp:patNum
    patJmp = min(patJmp, patNum - n);
    [minest, idx] = findMin(patWgt, patJmp);
    if(minest > patSize/2)
        break;
    end
    
    % index and weight
    gnpIdx = zeros(kNN, patJmp);
    gnpDis = zeros(kNN, patJmp);
    curPat = pat(:, idx);
    parfor m = 1:patJmp
        curPat_m = curPat(:, m);
        [gnpIdx_m, gnpDis_m] = vl_kdtreequery(kdTree, pat, curPat_m, ...
            'NumNeighbors', kNN);
        gnpIdx_m = gnpIdx_m';
        
        gnpIdx(:, m) = gnpIdx_m;
        
        gnpDis_m = (25/patSize^2)*gnpDis_m;
        gnpDis(:, m) = exp(-gnpDis_m');
    end
    
    % observations
    gupObv = cell(1, patJmp);
    gayPat = cell(1, patJmp);
    wamSat = cell(1, patJmp);
    for m = 1:patJmp
        gnpIdx_m = gnpIdx(:, m);
        
        gupObv_m = [patR(:,gnpIdx_m), patG(:,gnpIdx_m), patB(:,gnpIdx_m)];
        gayPat_m = 3*double(pat(1:patSize^2, gnpIdx_m));
        wamSat_m = [resPatR(:,gnpIdx_m), resPatG(:,gnpIdx_m), resPatB(:,gnpIdx_m)];
        
        gupObv{m} = gupObv_m;        
        gayPat{m} = gayPat_m;
        wamSat{m} = wamSat_m;
    end
    
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
        wamSat_m = wamSat{m};

        [newPat{m}, iter_m] = optADMM( gayPat_m, gupObv_m, ...
            Omega_m, mu, wamSat_m, para );
        iter(m) = length(iter_m);
    end
    iter = mean(iter);
    
    % weighted update results
    for m = 1:patJmp
        gnpIdx_m = gnpIdx(:, m);
        gnpDis_m = gnpDis(:, m);
        newPat_m = newPat{m};
        
        % update recover patches
        N = length(gnpIdx_m);
        resPatR(:, gnpIdx_m) = updateGupPat(newPat_m(:,1:N), ...
            resPatR(:, gnpIdx_m), gnpDis_m', patWgt(gnpIdx_m));
        resPatG(:, gnpIdx_m) = updateGupPat(newPat_m(:,N+1:2*N), ...
            resPatG(:, gnpIdx_m), gnpDis_m', patWgt(gnpIdx_m));
        resPatB(:, gnpIdx_m) = updateGupPat(newPat_m(:,2*N + 1:end), ...
            resPatB(:, gnpIdx_m), gnpDis_m', patWgt(gnpIdx_m));

        patWgt(gnpIdx_m) = patWgt(gnpIdx_m) + gnpDis_m';
    end
    
    % check convergence
    fprintf('%d out of %d, avg loops in opt %2.2f, min wgt: %2.2f \n', ...
        n, patNum, iter, minest);
end

resImg = patch2colorIm(resPatR, resPatG, resPatB, patIdx, patSize, imSize);

end

%% --------------------------------------------------------------
function [img] = patch2colorIm(patR, patG, patB, idx, patSize, imSize)

resImgR = patch2im(patR, idx, patSize, imSize);
resImgG = patch2im(patG, idx, patSize, imSize);
resImgB = patch2im(patB, idx, patSize, imSize);

img = cat(3, resImgR, resImgG, resImgB);

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

%% --------------------------------------------------------------
function [minest, idx] = findMin(wgt, num)
[~, idx] = sort(wgt, 'ascend');
idx = idx(1:num);

minest = wgt(idx(1));
end

