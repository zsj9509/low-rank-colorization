function [ L, obj ] = optADMM( W, O, Omega, mu, L, para )
% W: [m n] gray image
% O: [m 3n] observed positions
% Omega: [m 3n] observed values

% setup start and stop paras
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
else
    maxIter = 1000;
    tol = 1e-6;
end

% algorithm parameters
X = L;
Q = zeros(size(O));
lambda = min(20*numel(Omega)/nnz(Omega), 1000);
rho = 1e-3;
rho_max = 1e+10;

% temp matrix
[~, n] = size(W);
T = eye(n, n);
T = [T, T, T];
T = sparse(T');
A = sparse(T*T');
Omega = sparse(Omega);
tempC = W*T' + lambda*(Omega.*O);

% start loop
obj = zeros(maxIter, 1);
for i = 1:maxIter
    C = tempC + rho*X - Q;
    [ L, iter ] = conjGrad_admm( A, Omega, lambda, rho, C, L(:) );
    iter = length(iter);

    [U, S, V] = svt(L + Q/rho, mu/rho);
    X = U*S*V';

    Q = Q + rho*(L - X);

    obj(i) = (1/2)*sumsqr(L*T - W);
    obj(i) = obj(i) + (lambda/2)*sumsqr(Omega.*(L - O));
    obj(i) = obj(i) + mu*sum(diag(S));

    if(para.pnt == 1)
        fprintf('iter %d, obj %d, inner %d, rank %d \n', i, obj(i), iter, nnz(S));
    end
    
    if(rho < rho_max)
        rho = rho*1.25;
    end

    if(i > 3 && abs(obj(i) - obj(i - 1)) < tol)
        break;
    end
end

obj = obj(1:i);

end

