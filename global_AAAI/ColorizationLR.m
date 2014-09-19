function [L, S, X, obj] = ColorizationLR(Data, lambda, eta, tol, maxIter)

% observations
% indices
Omega = Data.Omega; 
% The partially observed R, G, B matrix
D = Data.D;
% the black&white matrix
B = 3 * Data.B;

[~, n] = size(B);
I = eye(n, n);
T = [I; I; I];
TT = T*T';

[m, n] = size(D);
I = eye(n);

if nargin < 2
    lambda = 1 / sqrt(m);
end

if nargin < 4
    tol = 1e-6;
elseif tol == -1
    tol = 1e-6;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
Y1 = zeros(size(D));
Y2 = Y1;

L = zeros( m, n);
X = zeros( m, n);

mu1 = 1; % this one can be tuned
mu2 = 1;

obj = zeros(maxIter, 1);
for iter = 1:maxIter        
    Zs = Y1/mu1 + D - L;
    S = ST_Shrinkage(Zs, Omega*lambda/mu1);
    
    Zl = ( Y1 + Y2 + mu1*(D - S) + mu2*X ) / (mu1 + mu2);
    [L,  sv] = SVShrinkage(Zl, 1/(mu1+mu2));
    
    X = (eta*TT + mu2*I) \ (eta*T*B' - Y2' + mu2*L');
    X = X';
    
    Z1 = (D - L - S);
    Z2 = (X - L);
    
    obj_iter = sum(sv(:));
    obj_iter = obj_iter + lambda*sum(abs(S(:) .* Omega(:)));
    obj_iter = obj_iter + (eta/2)*sumsqr(X*T - B);
    obj_iter = obj_iter + trace(Y1'*Z1);
    obj_iter = obj_iter + (mu1/2)*sumsqr(Z1);
    obj_iter = obj_iter + trace(Y2'*Z2);
    obj_iter = obj_iter + (mu2/2)*sumsqr(Z2);
    obj(iter) = obj_iter;
    
    fprintf('iter %d, obj: %d \n', iter,  obj_iter);
    
    Y1 = Y1 + mu1 * Z1;
    Y2 = Y2 + mu2 * Z2;
    
    if(iter > 3 && abs(obj(iter) - obj(iter - 1)) < tol)
        break;
    end
end

obj = obj(1:iter);
