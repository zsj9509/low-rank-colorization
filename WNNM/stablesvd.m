function [ U, S, V ] = stablesvd( X )

[m, n] = size(X);
if(m < n)
X = X';
end

try
    [U, S, V] = svd(X, 'econ');
catch
    save('err-svd.mat');
    [U,S,V] = lansvd(double(X),min(m, n),'L');
    U = single(U);
    V = single(V);
    S = single(S);
end

if(m < n)
    temp = U;
    U = V;
    V = temp;
end

end

