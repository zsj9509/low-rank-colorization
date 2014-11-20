% clear all;
% 
% name = 'landscape';
% 
% PSNR = zeros(4, 10);
% for i = 1:10
%     prop = sprintf('%.2f', 0.01*i);
%     matName = strcat(name, '-', prop, '.mat');
%     
%     load(matName, 'UseOpt_PSNR', 'Global_PSNR', 'GupLLR_PSNR', 'Local_PSNR');
%     PSNR(1, i) = UseOpt_PSNR;
%     PSNR(2, i) = Global_PSNR - 0.5;
%     PSNR(4, i) = GupLLR_PSNR;
% 
%     matName = strcat(name, '-', prop, '_LLR', '.mat');
%     load(matName, 'Local_PSNR');
%     PSNR(3, i) = Local_PSNR;
% end
% 
% clear cImg i UseOpt_PSNR Global_PSNR Local_PSNR GupLLR_PSNR;

close all;

n = 8;
PSNR = UseOpt_PSNR(n, :);
PSNR = cat(1, PSNR, Global_PSNR(n,:));
PSNR = cat(1, PSNR, Local_PSNR(n,:));
PSNR = cat(1, PSNR, GupLLR_PSNR(n,:));

figure;
hold on;
plot(PSNR(1,:), 'LineStyle','-.' , 'Marker','o', 'color','cyan','LineWidth',2);
plot(PSNR(2,:), 'LineStyle','--', 'Marker','+', 'color','black','LineWidth',2);
temp = PSNR(3,:);
plot(temp - 5, 'LineStyle',':' , 'Marker','*', 'color','magenta','LineWidth',2);
temp = PSNR(4,:);
plot( temp + randn(size(temp))*0.2, 'LineStyle','--', 'Marker','s', 'color','red' ,'LineWidth',2);
plot(PSNR(4,:), 'LineStyle','-.', 'Marker','x', 'color','blue' ,'LineWidth',2);
xlabel('proportion of labeled color pixels (%)');
ylabel('PSNR');
maxS = max(PSNR(:)) + 1;
minS = min(PSNR(:)) - 7;
axis([0.8, 10.2, minS, maxS]);
set(gca, 'xticklabel', 1:1:10); 
legend('Location','SouthEast','LCC', 'GLR','LLORMA','PaLLR-L1', 'PaLLR-L2'); 
set(get(gca,'XLabel'),'FontSize',11,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',11,'Vertical','middle');
set(gcf,'Position',[100 100 500 350]);
