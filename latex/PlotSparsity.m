surf(PSNR2);
xlabel('lambda');
ylabel('mu');
zlabel('PSNR');

set(gca, 'xticklabel', 0.01*2.^[0:4:12]);
set(gca, 'yticklabel', 0.01*2.^[0:2:15]);