function [ L, output ] = optADMM( W, O, Omega, mu, para )
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
L = repmat(W, 1, 3);
X0 = L;
X1 = X0;
alpha0 = 1;
alpha1 = 1;
Q = zeros(size(O));
if(isfield(para, 'lambda'))
    lambda = para.lambda*numel(Omega)/nnz(Omega);
else
    lambda = min(10*numel(Omega)/nnz(Omega), 400);
end
rho = 1e-2;
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
    hX = X1 + ((alpha0 - 1)/alpha1)*(X1 - X0);
    
    C = tempC + rho*hX - Q;
    [ L, iter ] = conjGrad_admm( A, Omega, lambda, rho, C, L(:) );
    iter = length(iter);

    [U, S, V] = svt(L + Q/rho, mu/rho);
    hX = U*S*V';

    Q = Q + rho*(L - hX);
    
    % acceleration updating
    X0 = X1;
    X1 = hX;
    
    halpha = (1 + sqrt(1 + alpha1^2))/2;
    alpha0 = alpha1;
    alpha1 = halpha;

    % check object value
    obj(i) = (1/2)*sumsqr(L*T - W);
    obj(i) = obj(i) + (lambda/2)*sumsqr(Omega.*(L - O));
    obj(i) = obj(i) + mu*sum(diag(S));

    if(para.pnt == 1)
        fprintf('iter %d, obj %d, inner %d, rank %d \n', i, obj(i), iter, nnz(S));
    end
    
    if(i > 3 && abs(obj(i) - obj(i - 1)) < tol)
        break;
    end
    
    % approximation
    if(rho < rho_max)
        rho = rho*1.15;
    end
end

output.obj = obj(1:i);
output.rank = nnz(S);

end

