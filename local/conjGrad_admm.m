function [ X, err ] = conjGrad_admm( A, Omega, lambda, rho, C, x, tol)

[m, n] = size(C);

b = C(:);
Omega = lambda*Omega;

if(~exist('x', 'var'))
    x = zeros(m*n, 1);
end

maxIter = ceil(sqrt(max(m, n)));
err = zeros(maxIter, 1);

% conjugate descent
r = fastAp( A, Omega, rho, x) - b;
p = - r;
for i = 1:maxIter
    Ap = fastAp( A, Omega, rho, p);
    alpha = (r'*r)/(p'*Ap);
    x = x + alpha*p;
    rk = r + alpha*Ap;
    beta = (rk'*rk)/(r'*r);
    p = - rk + beta*p;
    
    err(i) = sumsqr(rk);
    if(err(i) < tol)
        break;
    end
    r = rk;
end
err = err(1:i);

X = reshape(x, m, n);

end

%% --------------------------------------------------------------
function [Ap] = fastAp( A, Omega, rho, p)

[m, n] = size(Omega);
P = reshape(p, m, n);
Ap = P*A + (Omega .* P) + rho*P;

Ap = reshape(Ap, m*n, 1);

end

