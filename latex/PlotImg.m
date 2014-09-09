close all;

% figure;
% imshow(cImg);
% title('origin');

figure;
imshow(recover.Global);
title('global');

figure;
imshow(recover.UseOpt);
title('LCP');

figure;
imshow(recover.LLR);
title('LLR');

figure;
imshow(recover.GupLLR);
title('GupLLR');