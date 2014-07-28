function [ L, obj ] = optReweight( W, O, Omega, mu )
% W: [m n] gray image
% O: [m 3n] observed positions
% Omega: [m 3n] observed values

[~, n] = size(W);

L = repmat(W, 1, 3);
T = eye(n, n);
T = [T, T, T]';
T = sparse(T);
lambda = min(20*numel(Omega)/nnz(Omega), 1000);

Omega = sparse(Omega);
A = T*T';
C = W*T' + lambda*(Omega.*O);
LT = L*T;

tol = 1e-4;
maxIter = 500;
obj = zeros(1, maxIter);
for i = 1:maxIter
    [Q, invQ] = sqrtMatrix(L*L', 2e-3);
    
    B = mu*invQ;
    
    [L, iter] = conjGrad_rw(A, B, C, Omega, lambda, L(:));
    iter = length(iter);
    
    obj(i) = 0.5*sumsqr( LT - W );
    obj(i) = obj(i) + (lambda/2)*sumsqr(Omega.*(L - O));
    obj(i) = obj(i) + mu*trace(Q);
   
    fprintf('iter: %d obj %d, inner loop %d \n', i, obj(i), iter);
    if(i > 2 && abs(obj(i) - obj(i - 1)) < tol)
        break;
    end
   
end

obj = obj(1:i);

end

%% --------------------------------------------------------------
function [M, invM] = sqrtMatrix(X, epsilon)

X = X + epsilon*eye(size(X));
[E, D] = eig(X);
D = diag(D);

D = sqrt(D);
M = E*diag(D)*E';

D = 1./D;
invM = E*diag(D)*E';

end

