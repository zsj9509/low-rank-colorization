function [ L, resRank ] = optADMM( W, O, Omega, mu )

[m, n] = size(Omega);

T = eye(n/3, n/3);
T = [T, T, T];
T = T';

tol = 1e-8;

rho = 10;
R = zeros(m, n);
Q = zeros(m, n);
L = O;

WTt = W*T';
invT = T*T';
invT = invT + rho*eye(size(invT));
invT = inv(invT);
invT = sparse(invT);

maxIter = 1000;
err = zeros(maxIter, 1);
for i = 1:maxIter
    % update X
    X = min_X(WTt, rho*L, Q, R.*Omega, invT);
    % update L
    [U, S, V] = min_L(X, Q, mu, rho);
    L = U*S*V';
    
    % update dual R & S
    R = R + rho*(X - O).*Omega;
    Q = Q + rho*(X - L);
    
    err(i) = norm((X - O).*Omega) + norm(X - L);
    if(err(i) < tol)
        break;
    end
    
    fprintf('iter %d, error: %d, rank %d \n', i, err(i), nnz(S));        
end

resRank = nnz(S);

end

%% --------------------------------------------------------------
function X = min_X( WTt, rhoL, S, R_O, invT)

A = WTt + rhoL - S - R_O;

X = A*invT;

end

%% --------------------------------------------------------------
function [U, S, V] = min_L(X, S, mu, rho)

[U, S, V] = svd(X + S/rho, 'econ');

S = diag(S);
S = S - mu/rho;
S = S(S > 0);

U = U(:,1:length(S));
V = V(:,1:length(S));
S = diag(S);

end

