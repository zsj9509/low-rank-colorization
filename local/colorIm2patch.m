function [rPat, gPat, bPat] = colorIm2patch(cImg, patSize, sliding)

R = cImg(:,:,1);
G = cImg(:,:,2);
B = cImg(:,:,3);

rPat  = im2patch(R, patSize, sliding);
gPat  = im2patch(G, patSize, sliding);
bPat  = im2patch(B, patSize, sliding);

end