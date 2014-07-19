function [ L ] = optActive( W, O, Omega, mu )

[m, n] = size(W);
lambda = 10;
rho = initialRho(Omega, lambda);

T = eye(n, n);
T = [T, T, T]';

U = eye(m, min(m, 3*n));
V = eye(3*n, min(m, 3*n));
S = eye(min(m, 3*n));
L = U*S*V';

maxIter = 1000;
err = zeros(maxIter, 1);
for i = 1:maxIter    
    grad = (L*T - W)*T' + lambda*(L - O).*Omega;
    [U_G, S_G, V_G] = svd(L - grad, 'econ');
    S_G = sum(S_G(:) > mu);
    U_G = U_G(:, 1:S_G);
    V_G = V_G(:, 1:S_G);
    
    [U_A, ~] = qr([U_G, U], 0);
    [V_A, ~] = qr([V_G, V], 0);
    U_A = U_A(:, 1:min(size(U_A,2), size(V_A,2)));
    V_A = V_A(:, 1:min(size(U_A,2), size(V_A,2)));
    
    [U_S, S, V_S, inIter] = min_S(U_A, V_A, W, T, O, Omega, mu, lambda, rho);
    
    U = U_A*U_S;
    V = V_A*V_S;
    L = U*S*V';
    
    err(i) = objfunction(L*T - W, Omega.*(L - O), S, lambda, mu);
    
    fprintf('out %d \n', err(i));
end

end

%% --------------------------------------------------------------
function [U, S, V, err] = min_S(U_A, V_A, W, T, O, Omega, mu, lambda, rho)

sL = eye(size(U_A, 2));

maxIter = 1000;
err = zeros(maxIter, 1);
for i = 1:maxIter
    LT_W = U_A*sL*V_A'*T - W;
    L_O  = Omega.*(U_A*sL*V_A' - O);
    gradS = U_A'*LT_W*T'*V_A + lambda*U_A'*L_O*V_A;
    
    [U, S, V] = svt(sL - gradS/rho, mu/rho);
    sL = U*S*V';
    
    err(i) = objfunction(LT_W, L_O, S, lambda, mu);
    
    if(i > 2 && abs(err(i) - err(i - 1)) < 1e-6)
        break;
    end
    
    %fprintf('inner %d \n', err(i));
end
err = err(1:i);

end

%% --------------------------------------------------------------
function [obj] = objfunction(LT_W, L_O, S, lambda, mu)

obj = 0.5*sumsqr(LT_W);
obj = obj + lambda/2*sumsqr(L_O);
obj = obj + mu*sum(diag(S));

end

%% --------------------------------------------------------------
function rho = initialRho(Omega, lambda)
[~, pos] = max(sum(Omega, 1));
I = lambda*diag(Omega(:, pos));
I = I + 3*eye(size(I));
rho = svds(I, 1);
end

