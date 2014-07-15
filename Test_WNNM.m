clear all;

clnImg = imread('house.png');
% clnImg = rgb2gray(clnImg);
clnImg = double(clnImg)/255;

sigma = 10/255;

nosImg = clnImg + sigma*randn(size(clnImg));

% [Par]=ParSet(sigma);
% resImg =  WNNM_DeNoising( nosImg, clnImg, Par );
% 
% PSNR = 20*log10(1/sqrt(mean((resImg(:) - clnImg(:)).^2)));

[ resImg, PSNR ] = Denoising( nosImg, clnImg, 0.25 );