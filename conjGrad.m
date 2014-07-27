function [X, err] = conjGrad(A, B, C, Omega, lambda)
%% --------------------------------------------------------------

m = size(B, 2);
n = size(A, 1);

Omega = sparse(Omega(:));
Omega = diag(Omega);

L = kron(A', eye(m, m));
L = L + kron(eye(n, n), B);
L = L + lambda*Omega;

b = C(:);

if(~exist('x', 'var'))
    x = zeros(m*n, 1);
end

maxIter = m*n + min(m, n);
err = zeros(maxIter, 1);

clear A B C Omega lambda;

% conjugate descent
r = L*x - b;
p = - r;
for i = 1:maxIter
    alpha = (r'*r)/(p'*L*p);
    x = x + alpha*p;
    rk = r + alpha*L*p;
    beta = (rk'*rk)/(r'*r);
    p = - rk + beta*p;
    
    err(i) = norm(r, 2);
    if(err(i) < 1e-5)
        break;
    end
    r = rk;
end
err = err(1:i);

X = reshape(x, m, n);

end
