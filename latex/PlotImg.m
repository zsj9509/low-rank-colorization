close all;

% figure;
% imshow(cImg);
% title('origin');
i = 8;
t = 10;

figure;
imshow(recover{i,t}.Global);
title('GLR');

figure;
imshow(recover{i,t}.UseOpt);
title('LCC');

figure;
imshow(recover{i,t}.LLR);
title('LLR');

figure;
imshow(recover{i,t}.GupLLR);
title('GLLR');