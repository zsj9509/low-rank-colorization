close all;

figure;
hold on;
p = plot(SSIM, 'LineStyle','-' , 'Marker','o' ,'LineWidth',2);
axis([0.8, 9.2, 0.98, 0.988]);
set(gca, 'xticklabel', 6:2:22); 

ylabel('SSIM');
xlabel('Patch Size');
