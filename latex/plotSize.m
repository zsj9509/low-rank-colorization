%% ---------------------------------------------------------------
close all;

clear all;

load('castlePatSize.mat', 'PSNR');
patSize(1,:) = PSNR(1:9);
clear PSNR;
% load('landscapePatSize.mat');
% patSize(2,:) = PSNR(1:9);
% clear PSNR;

figure;
hold on;
p = plot(patSize(1,:), 'LineStyle','--' , 'Marker','o' ,'LineWidth',2);
% p = plot(patSize(2,:), 'LineStyle','-' , 'Marker','+' ,'LineWidth',2,'color', 'red');
axis([0.8, 9.2, 33, 43]);
set(gca, 'xticklabel', 4:2:20); 

ylabel('PSNR');
xlabel('Patch Size');
set(get(gca,'XLabel'),'FontSize',11,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',11,'Vertical','middle');

%% ---------------------------------------------------------------
% clear all;
% load('castleGupSize.mat', 'PSNR');
% patSize(1,:) = PSNR(1:11);
% clear PSNR;
% % load('landscapeGupSize.mat', 'PSNR');
% % patSize(2,:) = PSNR(1:11);
% % clear PSNR;
% 
% figure;
% hold on;
% p = plot(patSize(1,:), 'LineStyle','--' , 'Marker','o' ,'LineWidth',2);
% % p = plot(patSize(2,:), 'LineStyle','-' , 'Marker','+' ,'LineWidth',2, 'color', 'red');
% axis([0.8, 11.2, 39, 42]);
% set(gca, 'xticklabel', 10:10:110); 
% 
% ylabel('PSNR');
% xlabel('Group Size');
% 
% ylabel('PSNR');
% xlabel('Group Size');
% set(get(gca,'XLabel'),'FontSize',11,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',11,'Vertical','middle');