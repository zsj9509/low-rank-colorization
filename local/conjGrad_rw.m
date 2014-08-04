function [X, err] = conjGrad_rw(A, B, C, Omega, lambda, x)
%% --------------------------------------------------------------

m = size(B, 2);
n = size(A, 1);
b = C(:);

if(~exist('x', 'var'))
    x = zeros(m*n, 1);
end

maxIter = ceil(sqrt(min(m, n)));
err = zeros(maxIter, 1);

% conjugate descent
r = fastLp(A, B, Omega, lambda, x) - b;
p = - r;
for i = 1:maxIter
    Lp = fastLp(A, B, Omega, lambda, p);
    alpha = (r'*r)/(p'*Lp);
    x = x + alpha*p;
    rk = r + alpha*Lp;
    beta = (rk'*rk)/(r'*r);
    p = - rk + beta*p;
    
    err(i) = norm(rk, 2);
    if(err(i) < 1e-6)
        break;
    end
    r = rk;
end
err = err(1:i);

X = reshape(x, m, n);

end

%% ---------------------------------------------------------------
function [x] = fastLp(A, B, Omega, lambda, x)

m = size(B, 2);
n = size(A, 1);

X = reshape(x, m, n);
X1 = X*A;
X2 = B*X;
X3 =lambda*(Omega.*X);
X = X1 + X2 + X3;

x = X(:);

end