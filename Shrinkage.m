function [Z] = Shrinkage(Y, tau)

Z = max(Y - tau, 0) + min(Y + tau, 0);

end