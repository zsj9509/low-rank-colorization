function [L, S, X, iter] = ColorizationLR(Data, lambda, eta, tol, maxIter)

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
    tol = 1e-8;
elseif tol == -1
    tol = 1e-8;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
Y1 = [B B B] + D;
norm_two = lansvd(Y1, 1, 'L');
norm_inf = norm( Y1(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y1 = Y1 / dual_norm;

Y2 = Y1;

S = zeros( m, n);

L = zeros( m, n);
X = zeros( m, n);


mu1 = 1.25/norm_two; % this one can be tuned
mu2 = 1.25/norm_two;
mu_bar = mu1 * 1e15;
rho = 1.5;          % this one can be tuned
%rho = 3;

iter = 0;
total_svd = 0;
converged = false;
sv = 10;

while ~converged       
    iter = iter + 1;
    
    Zs = Y1/mu1 + D - L;
    S = ST_Shrinkage(Zs, Omega*lambda/mu1);
    
    Zl = ( Y1 + Y2 + mu1*(D - S) + mu2*X ) / (mu1 + mu2);
    [L, sv] = SVShrinkage(Zl, 1/(mu1+mu2), sv, 0);
    
    X = (eta*TT + mu2*I) \ (eta*T*B' - Y2' + mu2*L');
    X = X';
    
    Z1 = (D - L - S);
    Y1 = Y1 + mu1 * Z1;
    Z2 = (X - L);
    Y2 = Y2 + mu2 * Z2;

    total_svd = total_svd + 1;
    
    mu1 = min(mu1 * rho, mu_bar);
    mu2 = min(mu2 * rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z1, 'fro') / norm(B, 'fro') + norm(Z2) / norm(B, 'fro');
    if stopCriterion < tol
        converged = true;
    end    
    
    if mod( total_svd, 10) == 0
        disp(['#svd ' num2str(total_svd) ' r(L) ' num2str(rank(L))...
            ' |S|_0 ' num2str(length(find(abs(S)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);
    end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end
