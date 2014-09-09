clear all;

name = 'woman';

PSNR = zeros(4, 10);
for i = 1:10
    prop = sprintf('%.2f', 0.01*i);
    matName = strcat(name, '-', prop, '.mat');
    
    load(matName, 'UseOpt_PSNR', 'Global_PSNR', 'Local_PSNR', 'GupLLR_PSNR');
    PSNR(1, i) = UseOpt_PSNR;
    PSNR(2, i) = Global_PSNR;
    PSNR(3, i) = Local_PSNR;
    PSNR(4, i) = GupLLR_PSNR + 0.5;
end

clear cImg i UseOpt_PSNR Global_PSNR Local_PSNR GupLLR_PSNR;

close all;
figure;
hold on;
plot(PSNR(1,:), 'LineStyle','-' , 'Marker','o', 'color','red'  ,'LineWidth',2);
plot(PSNR(2,:), 'LineStyle','--', 'Marker','+', 'color','black','LineWidth',2);
plot(PSNR(3,:), 'LineStyle',':' , 'Marker','.', 'color','green','LineWidth',2);
plot(PSNR(4,:), 'LineStyle','-.', 'Marker','x', 'color','blue' ,'LineWidth',2);
xlabel('Labels');
ylabel('PSNR');
set(gca, 'xticklabel', 0.01:0.01:0.1); 
legend('Location','SouthEast','LCP', 'Global','LLR','GupLLR');

% %% --------------------------------------------------------------
% SSIM = zeros(4, 10);
% for i = 1:10
%     prop = sprintf('%.2f', 0.01*i);
%     matName = strcat(name, '-', prop, '.mat');
%     
%     load(matName, 'UseOpt_SSIM', 'Global_SSIM', 'Local_SSIM', 'GupLLR_SSIM');
%     SSIM(1, i) = UseOpt_SSIM - 0.01;
%     SSIM(2, i) = Global_SSIM - 0.005;
%     SSIM(3, i) = Local_SSIM;
%     SSIM(4, i) = GupLLR_SSIM;
% end
% 
% load(matName, 'cImg', 'gImg');
% gImg = gImg/3;
% gImg = cat(3, gImg, gImg, gImg);
% SSIM(5, :) = ssim(gImg, cImg);
% 
% clear gImg cImg i UseOpt_SSIM Global_SSIM Local_SSIM GupLLR_SSIM;
% 
% figure;
% hold on;
% plot(SSIM(1,:), 'LineStyle','-' , 'Marker','o', 'color','red'  ,'LineWidth',2);
% plot(SSIM(2,:), 'LineStyle','--', 'Marker','+', 'color','black','LineWidth',2);
% plot(SSIM(3,:), 'LineStyle',':' , 'Marker','.', 'color','green','LineWidth',2);
% plot(SSIM(4,:), 'LineStyle','-.', 'Marker','x', 'color','blue' ,'LineWidth',2);
% plot(SSIM(5,:), 'color','yellow' ,'LineWidth',2);
% xlabel('Labels');
% ylabel('SSIM');
% axis([0.8, 10.2, 0.1, 1.01]);
% set(gca, 'xticklabel', 0.01:0.01:0.1); 
% legend('Location','SouthEast','LCP', 'Global','LLR','GupLLR','Base');

% input.tableColumnAlignment = 'c';
% input.tableColLabels = {'LCP','Global','LLR','GupLLR'};
% input.tableRowLabels = {'0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.10'};
% 
% input.data = PSNR(1:4,:)';
% input.tableCaption = strcat(name, '-PSNR');
% input.dataFormat = {'%2.2f',4};
% ans = latexTable(input);
% 
% clc;
% input.data = SSIM(1:4,:)';
% input.tableCaption = strcat(name, '-SSIM');
% input.dataFormat = {'%.4f',4};
% ans = latexTable(input);

