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
L = zeros(size(W)); L = repmat(L, 1, 3);
X0 =  zeros(size(L));
hX = X0;
alpha0 = 1;
alpha1 = 1;
Q0 = zeros(size(O));
hQ = Q0;
if(isfield(para, 'lambda'))
    lambda = para.lambda*numel(Omega)/nnz(Omega);
else
    lambda = min(10*numel(Omega)/nnz(Omega), 400);
end
rho = 1e-3;

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
    C = tempC + rho*hX - hQ;
    [ L, iter ] = conjGrad_admm( A, Omega, lambda, rho, C, L(:), max(1/i^3,1e-10) );
    iter = length(iter);

    [U, S, V] = svt(L + hQ/rho, mu/rho);
    X1 = U*S*V';

    tempLX1 = L - X1;
    % check object value
    obj(i) = (1/2)*sumsqr(L*T - W);
    obj(i) = obj(i) + (lambda/2)*sumsqr(Omega.*(L - O));
    obj(i) = obj(i) + mu*sum(diag(S));
    % obj(i) = obj(i) + trace(hQ'*tempLX1);
    % obj(i) = obj(i) + (rho/2)*sumsqr(tempLX1);

    if(para.pnt == 1)
        fprintf('iter %d, obj %d, inner %d, rank %d \n', i, obj(i), iter, nnz(S));
    end
    
    if(i > 3 && abs(obj(i) - obj(i - 1)) < tol)
        break;
    end
    
    Q1 = hQ + rho*tempLX1;
    % acceleration updating
    hX = X1 + ((alpha0 - 1)/alpha1)*(X1 - X0);
    X0 = X1;
    
    hQ = Q1 + ((alpha0 - 1)/alpha1)*(Q1 - Q0);
    Q0 = Q1;
    
    if(para.acc == 1)
        halpha = (1 + sqrt(1 + alpha1^2))/2;
        alpha0 = alpha1;
        alpha1 = halpha;
    else
        alpha0 = 1;
        alpha1 = 1;
    end
    
    rho = rho*1.25;
end

output.obj = obj(1:i);
output.rank = nnz(S);

end

