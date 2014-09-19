function [Z, S] = SVShrinkage(Y, tau)


[U, S, V] = svt(Y,  tau);
Z = U*S*V';

end