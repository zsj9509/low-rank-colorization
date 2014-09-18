close all;

% figure;
% imshow(cImg);
% title('origin');

figure;
imshow(recover.Global);
title('GLR');

figure;
imshow(recover.UseOpt);
title('LCC');

figure;
imshow(recover.LLR);
title('LLR');

figure;
imshow(recover.GupLLR);
title('GLLR');