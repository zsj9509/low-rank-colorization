DUseOpt = recover.UseOpt;
DUseOpt = DUseOpt - cImg;
DGlobal = recover.Global;
DGlobal = DGlobal - cImg;
DLLR = recover.LLR;
DLLR = DLLR - cImg;
DGupLLR = recover.GupLLR;
DGupLLR = DGupLLR - cImg;

DUseOpt = abs(DUseOpt);
DGlobal = abs(DGlobal);
DGupLLR = abs(DGupLLR);
DLLR = abs(DLLR);

DUseOpt = sum(DUseOpt, 3);
DGlobal = sum(DGlobal, 3);
DGupLLR = sum(DGupLLR, 3);
DLLR = sum(DLLR, 3);

normal = max([max(DUseOpt(:)), max(DGlobal(:)), max(DGupLLR(:))]);

DUseOpt = DUseOpt/normal;
DUseOpt = DUseOpt*255;
DUseOpt = uint8(DUseOpt);
DGlobal = DGlobal/normal;
DGlobal = DGlobal*255;
DGlobal = uint8(DGlobal);
DGupLLR = DGupLLR/normal;
DGupLLR = DGupLLR*255;
DGupLLR = uint8(DGupLLR);
DLLR = DLLR*255;
DLLR = uint8(DLLR);

close all;
figure;
imagesc(DUseOpt);
caxis([0 255]);
colorbar('location', 'westoutside'); 
axis off
figure;
imagesc(DGlobal);
caxis([0 255]);
colorbar('location', 'westoutside'); 
axis off
figure;
imagesc(DGupLLR);
caxis([0 255]);
colorbar('location', 'westoutside'); 
axis off
figure;
imagesc(DLLR);
caxis([0 255]);
colorbar('location', 'westoutside'); 
axis off
