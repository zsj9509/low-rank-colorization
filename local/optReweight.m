function [ L, obj ] = optReweight( W, O, Omega, mu, para )
% W: [m n] gray image
% O: [m 3n] observed positions
% Omega: [m 3n] observed values

% setup start and stop paras
[m, n] = size(W);
if(exist('para', 'var'))
    if(isfield(para, 'tol'))
        tol = para.tol;
    else
        tol = 1e-6;
    end
    
    if(isfield(para, 'maxIter'))
        maxIter = para.maxIter;
    else
        maxIter = 1000;
    end
    
    if(isfield(para, 'L'))
        L = para.L;
    else
        L = eye(m, 3*n);
    end
else
    L = eye(m, 3*n);
    maxIter = 1000;
    tol = 1e-6;
end

T = eye(n, n);
T = [T, T, T]';
T = sparse(T);
lambda = min(20*numel(Omega)/nnz(Omega), 1000);

Omega = sparse(Omega);
A = T*T';
C = W*T' + lambda*(Omega.*O);

obj = zeros(maxIter, 1);
for i = 1:maxIter
    [Q, invQ] = sqrtMatrix(L*L', 1e-3);
    
    B = mu*invQ;
    
    [L, iter] = conjGrad_rw(A, B, C, Omega, lambda, L(:));
    iter = length(iter);
    
    obj(i) = 0.5*sumsqr( L*T - W );
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

