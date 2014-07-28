function [ L1, obj ] = optProximal( W, O, Omega, mu, para )
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
    
    if(isfield(para, 'L'))
        L0 = para.L;
        L1 = repmat(W, 1, 3);
    else
        L0 = zeros(size(O));
        L1 = L0;
    end
else
    L0 = zeros(size(O));
    L1 = repmat(W, 1, 3);
    maxIter = 1000;
    tol = 1e-6;
end

t0 = 1;
t1 = 1;

% algorithm matrix
T = eye(size(W, 2));
T = [T, T, T]';
T = sparse(T);
TTt = T*T';
WTt = W*T';

lambda = min(20*numel(Omega)/nnz(Omega), 1000);
rho = initialRho(Omega, lambda);

obj    = zeros(maxIter, 1);
obj(1) = objectValue(L0, T, W, O, Omega, 0, lambda, mu);
for i = 2:maxIter       
    acc = 1;
    if(acc)
        Y = L1 + (t0 - 1)/t1 * (L1 - L0);
    else
        Y = L0;
    end
    
    gradY = Y*TTt - WTt + lambda*(Y - O).*Omega;
    [U, S, V] = svt(Y - gradY/rho, mu/rho);
    newL = U*S*V';
    L0 = L1;
    L1 = newL;

    obj(i) = objectValue(newL, T, W, O, Omega, S, lambda, mu);
    
    t0 = t1;
    t1 = 1 + sqrt(1 + 4*t0^2)/2;
    
    fprintf('iter %d, obj %d, rank %d, violate %d \n', i, obj(i), nnz(S), ...
        sumsqr(Omega.*(L0 - O)));
    
    if(i > 3 && abs(obj(i) - obj(i - 1)) < tol)
        break;
    end
end

obj = obj(2:i);

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

