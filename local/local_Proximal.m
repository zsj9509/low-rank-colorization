function [ L0, obj ] = local_Proximal( W, O, Omega, mu )

tol = 1e-3;
[m, n] = size(O);

L0 = zeros(m, n);
L1 = zeros(m, n);

t0 = 1;
t1 = 1;

T = eye(size(W, 2));
T = [T, T, T]';
T = sparse(T);

lambda = 5*numel(Omega)/nnz(Omega);
rho = initialRho(Omega, lambda);
maxIter = 1000;
obj    = zeros(maxIter, 1);
obj(1) = objectValue(L0, T, W, O, Omega, 0, lambda, mu);

for i = 2:maxIter    
%     gradL = (L0*T - W)*T' + lambda*(L0 - O).*Omega;
%     while(true)        
%         [U, newS, V] = proximal_L(L0, gradL, mu, rho);
%         newL = U*newS*V';
%         
%         new_obj = objectValue(newL, T, W, O, Omega, newS, lambda, mu);
%         pro_obj = trace(gradL'*(newL - L0)) + (0.5*rho)*norm(newL - L0, 'fro')^2;
%         
%         if(new_obj < obj(i - 1) + pro_obj)
%             L0 = newL;
%             obj(i) = new_obj;
%             break;
%         else
%             rho = rho*2;
%         end
%     end    
    
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

    obj(i) = objectValue(L0, T, W, O, Omega, S, lambda, mu);
    
    t0 = t1;
    t1 = 1 + sqrt(1 + 4*t0^2)/2;
    
    fprintf('iter %d, obj %d, rank %d, violate %d \n', i, obj(i), nnz(S), sumsqr(Omega.*(L0 - O)));
    
    if(i > 3 && abs(obj(i) - obj(i - 1)) < tol)
        break;
    end
end

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

