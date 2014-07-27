function [ rImg ] = colorUseOpt( gI, Obvs )



colorIm = Obvs(:,:,1);
colorIm = double(colorIm > 0);

idx = find(Obvs(:,:,1) > 0);
cI = zeros(size(Obvs));
temp1 = gI;
temp2 = Obvs(:,:,1);
% temp2 = temp2 + randn(size(temp2))*0.02;
temp1(idx) = temp2(idx);
cI(:,:,1) = temp1;

temp1 = gI;
temp2 = Obvs(:,:,2);
% temp2 = temp2 + randn(size(temp2))*0.02;
temp1(idx) = temp2(idx);
cI(:,:,2) = temp1;

temp1 = gI;
temp2 = Obvs(:,:,3);
% temp2 = temp2 + randn(size(temp2))*0.02;
temp1(idx) = temp2(idx);
cI(:,:,3) = temp1;

gI = cat(3, gI, gI, gI);
sgI=rgb2ntsc(gI);
scI=rgb2ntsc(cI);
   
ntscIm(:,:,1)=sgI(:,:,1);
ntscIm(:,:,2)=scI(:,:,2);
ntscIm(:,:,3)=scI(:,:,3);

rImg = getColorExact(colorIm,ntscIm);

end

