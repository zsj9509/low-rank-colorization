function [img] = patch2colorIm(patR, patG, patB, idx, patSize, imSize)

resImgR = patch2im(patR, idx, patSize, imSize);
resImgG = patch2im(patG, idx, patSize, imSize);
resImgB = patch2im(patB, idx, patSize, imSize);

img = cat(3, resImgR, resImgG, resImgB);

end