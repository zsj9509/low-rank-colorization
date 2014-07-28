function [ L1, rank, obj ] = optProximal( W, O, Omega, mu, initL )
% W: [m n] gray image
% O: [m 3n] observed positions
% Omega: [m 3n] observed values

[m, n] = size(O);

if(exist('initL','var'))
    L0 = initL;
    L1 = initL;
else
    L0 = zeros(m, n);
    L1 = zeros(m, n);
end

t0 = 1;
t1 = 1;

T = eye(size(W, 2));
T = [T, T, T]';
T = sparse(T);

lambda = min(20*numel(Omega)/nnz(Omega), 1000);
rho = initialRho(Omega, lambda);
maxIter = 10000;
obj    = zeros(maxIter, 1);
obj(1) = objectValue(L0, T, W, O, Omega, 0, lambda, mu);

tol = 1e-6;

for i = 2:maxIter       
    acc = 1;
    if(acc)
        Y = L1 + (t0 - 1)/t1 * (L1 - L0);
    else
        Y = L0;
    end
    
    gradY = (Y*T - W)*T' + lambda*(Y - O).*Omega;
    [U, S, V] = proximal_L(Y, gradY, mu, rho);
    newL = U*S*V';
    L0 = L1;
    L1 = newL;

    obj(i) = objectValue(newL, T, W, O, Omega, S, lambda, mu);
    
    t0 = t1;
    t1 = 1 + sqrt(1 + 4*t0^2)/2;
    
    fprintf('iter %d, obj %d, rank %d, violate %d \n', i, obj(i), nnz(S), sumsqr(Omega.*(L0 - O)));
    
    if(i > 3 && abs(obj(i) - obj(i - 1)) < tol)
        break;
    end
end

rank = nnz(S);

obj = obj(2:i);

end

%% --------------------------------------------------------------
function [U, S, V] = proximal_L(L, gradL, mu, rho)

Z = L - gradL/rho;
[U, S, V] = svd(Z, 'econ');

S = diag(S);
S = S - mu/rho;
S = S( S > 0 );

U = U(:, 1:length(S));
V = V(:, 1:length(S));
S = diag(S);

end

%% --------------------------------------------------------------
function obj = objectValue(L, T, W, O, Omega, S, lambda, mu)

obj = 0.5*sumsqr(L*T - W);
obj = obj + 0.5*lambda*sumsqr(Omega.*(L - O));
obj = obj + mu*sum(diag(S));

end

%% --------------------------------------------------------------
function rho = initialRho(Omega, lambda)
[~, pos] = max(sum(Omega, 1));
I = lambda*diag(Omega(:, pos));
I = I + 3*eye(size(I));
rho = svds(I, 1);
end

