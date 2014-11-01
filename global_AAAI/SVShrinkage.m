function [Z, sv] = SVShrinkage(Y, tau, sv, choose)
%Z = argmin_Z  tau ||Z||_* + 0.5 ||Z - Y||_F^2

[m n] = size(Y);

if nargin < 4
    choose = choosvd(n, sv);
end

% addpath PROPACK;


if choose == 1
    [U S V] = lansvd(Y, sv, 'L');
else
    [U S V] = svd(Y, 'econ');
end

diagS = diag(S);
svp = length(find(diagS > tau));

if svp < sv
    sv = min(svp + 1, n);
else
    sv = min(svp + round(0.05*n), n);
end

Z = U(:, 1:svp) * diag(diagS(1:svp) - tau) * V(:, 1:svp)';    

end