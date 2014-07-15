function [IMout, Weight] = patch2im(patches, idx, patchSize, imSize)

count = 1;
Weight=zeros(imSize);
IMout = zeros(imSize);
[rows,cols] = ind2sub(imSize-patchSize+1, idx);
for i  = 1:length(cols)
    col = cols(i); 
    row = rows(i);
    block = reshape(patches(:,count), [patchSize, patchSize]);
    IMout(row:row+patchSize-1,col:col+patchSize-1) = ...
        IMout(row:row+patchSize-1,col:col+patchSize-1)+block;
    Weight(row:row+patchSize-1,col:col+patchSize-1) = ...
        Weight(row:row+patchSize-1,col:col+patchSize-1)+ones(patchSize);
    count = count+1;
end;

end