function [Z] = IndexedShrinkage(Y, Omega, tau)
%Z = argmin_Z   tau * ||Z \circ Omega||_1 + 0.5 * ||Z - Y||_F^2

Z = max(Y - tau, 0) + min(Y + tau, 0);
Z(~Omega) = Y(~Omega);

end