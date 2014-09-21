function [L S X iter] = ColorizationLL(Data, lambda, eta, tolerance, r)

if nargin < 5
    [G] = LocalColorConsistency(Data.D, Data.B, Data.Omega, tolerance);
else
    [G] = LocalColorConsistency(Data.D, Data.B, Data.Omega, tolerance, r);
end

O = double(G > 0);
Data.Omega = O;
Data.D = G;
[L S X iter] = ColorizationLR(Data, lambda, eta);

end