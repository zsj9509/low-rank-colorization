function [patches, idx, pos] = im2patch(I, blkSize, slidingDis)

idxMat = zeros(size(I) - blkSize + 1);
% take blocks in distances of 'slidingDix', 
% but always take the first and last one (in each row and column).
idxMat([1:slidingDis:end-1,end],[1:slidingDis:end-1,end]) = 1; 
idx = find(idxMat);
[rows, cols] = ind2sub(size(idxMat),idx);
patches = zeros(blkSize^2,length(idx));
for i = 1:length(idx)
    currBlock = I(rows(i):rows(i) + blkSize -1,...
        cols(i):cols(i)+ blkSize -1);
    patches(:,i) = currBlock(:);
end

% pos(1):cols  pos(2):rows
pos = double(cat(1, cols', rows'));
