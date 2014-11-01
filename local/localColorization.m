function [ resImg ] = localColorization( gImg, Obvs, mu, patPara)

kNN = patPara.kNN;
patSize = patPara.patSize;
sliding = patPara.sliding;
epsilon = patPara.epsilon;
rho = patPara.rho;

% image to patch
imSize = size(gImg);
[pat, patIdx, patPos] = im2patch(gImg, patSize, sliding);
[patR, patG, patB] = colorIm2patch(Obvs, patSize, sliding);

patNum = length(patIdx);
if(isfield(patPara, 'rand') && patPara.rand == 1)
    idx = randperm(patNum);
    pat = pat(:,idx); 
    patIdx = patIdx(idx); 
    patPos = patPos(:, idx);
    patR = patR(:, idx);
    patG = patG(:, idx);
    patB = patB(:, idx);
    clear idx;
end

% build up kd tree
pat = augPatch(pat/3, rho*patPos, imSize, 2);
kdTree = vl_kdtreebuild(pat);

resPatR = zeros(size(patR));
resPatG = zeros(size(patG));
resPatB = zeros(size(patB));

patWgt = zeros(1, size(pat, 2));
patJmp = 400;
for n = 1:patJmp:patNum
    patJmp = min(patJmp, patNum - n);
    [minest, idx] = findMin(patWgt, patJmp);
    if(minest > epsilon)
        break;
    end
    
    % index and weight
    gnpIdx = zeros(kNN, patJmp);
    gnpDis = zeros(kNN, patJmp);
    curPat = pat(:, idx);
    parfor m = 1:patJmp
    % for m = 1:patJmp
        curPat_m = curPat(:, m);
        [gnpIdx_m, gnpDis_m] = vl_kdtreequery(kdTree, pat, curPat_m, ...
            'NumNeighbors', kNN);
        gnpIdx_m = gnpIdx_m';
        
        gnpIdx(:, m) = gnpIdx_m;
        
        gnpDis_m = (25/patSize^2)*gnpDis_m;
        gnpDis_m = exp(-gnpDis_m');
        gnpDis(:, m) = gnpDis_m;
    end
    
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
    
    % optimization
    para.maxIter = 200;
    para.tol = 1e-6;    
    para.pnt = 0;
    para.acc = 1;
    if(isfield(patPara, 'lambda'))
        para.lambda = patPara.lambda;
    end
    
    newPat = cell (1, patJmp);
    iter   = zeros(1, patJmp); 
    rank   = zeros(1, patJmp);
    parfor m = 1:patJmp
    % for m = 1:patJmp
        Omega_m = double(gupObv{m} ~= -1);
        Omega_m = sparse(Omega_m);
        gayPat_m = gayPat{m};
        gupObv_m = gupObv{m};
        
        if(nnz(gupObv_m) < 2)
            continue;
        end

        [newPat_m, out] = optADMM( gayPat_m, gupObv_m, Omega_m, mu, para );
        
        newPat_m(newPat_m < 0) = 0;
        newPat_m(newPat_m > 1) = 1;
        
        iter(m) = length(out.obj);
        rank(m) = out.rank;
        newPat{m} = newPat_m;
    end
    iter = mean(iter);
    rank = mean(rank);
    
    % weighted update results
    for m = 1:patJmp
        gupObv_m = gupObv{m};
        if(nnz(gupObv_m) < 2)
            continue;
        end
        
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
    if(patPara.pnt == 1)
        fprintf('%d out of %d, avg loops %2.2f rank %2.2f, min wgt: %2.2f \n', ...
        n, patNum, iter, rank, minest);
        fprintf('nnz: %1.2f\n', nnz(patWgt)/numel(patWgt));
    end
    
    if(isfield(patPara, 'cImg'))
        resImg = patch2colorIm(resPatR, resPatG, resPatB, patIdx, patSize, imSize);
        PSNR = psnr(resImg, patPara.cImg);
        fprintf('PSNR %2.4f , nnz: %d\n', PSNR, nnz(patWgt)/numel(patWgt));
    end
end

resImg = patch2colorIm(resPatR, resPatG, resPatB, patIdx, patSize, imSize);

end

%% ------------------------------------------------------------------------
function [pat] = updateGupPat(newPat, oldPat, newWgt, oldWgt)

dim = size(newPat, 1);

newWgt = repmat(newWgt, dim, 1);
oldWgt = repmat(oldWgt, dim, 1);

newPat = newPat.*newWgt;
oldPat = oldWgt.*oldPat;

pat = newPat + oldPat;
pat = pat ./ (newWgt + oldWgt);

end

%% ------------------------------------------------------------------------
function [minest, idx] = findMin(wgt, num)
[~, idx] = sort(wgt, 'ascend');
idx = idx(1:num);

minest = wgt(idx(1));
end

