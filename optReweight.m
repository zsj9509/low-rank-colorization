function [ L ] = optReweight( W, O, Omega, mu )

[~, n] = size(W);

L = repmat(W, 1, 3);
T = eye(n, n);
T = [T, T, T]';
T = sparse(T);
lambda = min(20*numel(Omega)/nnz(Omega), 1000);

maxIter = 1000;
for i = 1:maxIter
    [Q, invQ] = sqrtMatrix(L'*L, 1e-4);
    
    A = T*T';
    B = mu*invQ;
    C = W*T' + lambda*(Omega.*O);
    
    [L, iter] = conjGrad(A, B, C, Omega, lambda);
    iter = length(iter);
    
    obj = 0.5*norm(L*T - W)^2;
    obj = obj + (lambda/2)*norm(Omega.*(L - O));
    obj = obj + mu/2*(trace(Q) + trace(L'*invQ*L));
    
    fprintf('iter: %d obj %d, inner loop %d', i, obj, iter);
end

end

%% --------------------------------------------------------------
function [M, invM] = sqrtMatrix(X, epsilon)

I = eye(size(X));
X = X + epsilon*I;
[E, D] = eig(X);
D = diag(D);

D = sqrt(D);
M = E*diag(D)*E';

D = 1./D;
invM = E*diag(D)*E';

end

