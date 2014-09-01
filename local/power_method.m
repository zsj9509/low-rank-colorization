function [ U, S, V, err ] = power_method( A, R, maxIter)

if(~exist('tol','var'))
    tol = 1e-4;
end

Y = A*R;
[Q, ~] = qr(Y, 0);
err = zeros(maxIter, 1);
for i = 1:maxIter
    Y = A*(A'*Q);
    % Y = AAt*Q;
    [iQ, ~] = qr(Y, 0);
    err_i = norm(Q(1,:) - iQ(1,:), 'fro');
    err(i) = err_i;
    if(err_i < tol)
        break;
    end
    
    % fprintf('iter %d, err: %f\n', i, err_i);
    Q = iQ;
end
err = err(1:i);

B = Q'*A;
[mB, nB] = size(B);
if(mB > nB)
    [hU, S, V] = svd(B, 'econ');
else
    [hU, S, V] = svd(B', 'econ');
    temp = hU;
    hU = V;
    V = temp;
end
U = Q*hU;

end

