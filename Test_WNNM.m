clear all;

clnImg = imread('hotel.jpg');
clnImg = rgb2gray(clnImg);
clnImg = double(clnImg);

sigma = 0.1;

nosImg = clnImg + sigma*randn(size(clnImg));

[Par]=ParSet(sigma);
[E_Img, PSNR]   =  WNNM_DeNoising( nosImg, clnImg, Par );