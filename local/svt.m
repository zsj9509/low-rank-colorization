function [U, S, V] = svt(Z, mu)

[U, S, V] = stableSvd(Z);

S = diag(S);
S = S - mu;
S = S( S > 0 );

U = U(:, 1:length(S));
V = V(:, 1:length(S));
S = diag(S);

end

%% --------------------------------------------------------------
function [U, S, V] = stableSvd(Z)

[m, n] = size(Z);

if(m < n)
    Z = Z';
end

try
    [U, S, V] = svd(Z, 'econ');
catch
    [U, S, V] = lansvd(Z, min([m,n]), 'L');
end

if(m < n)
    temp = U;
    U = V;
    V = temp;
end

end