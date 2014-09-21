function [G] = LocalColorConsistency(D, B, Omega, tau, r)


if nargin < 4
    tau = 0.05;
end
if nargin < 5
    Omega = (D > 0);
    [size1, size2] = size(Omega);
    observed = sum(sum(Omega)) / size1 / size2;
    radius = sqrt(2 / observed);
    r = ceil(radius);
    if observed > 0.1
        r = 5;
    end
    if observed < 0.005
        r = 20;
    end
end

[M N] = size(B);
G = zeros(M, 3*N);

% we want to label each unlabeled (i, j)
for i = 1: M
    for j = 1: N
        % if (i, j) is already labeled, then skip it
        if Omega(i, j) == 1
            continue;  
        end
        isFindSimilarNeighbor = 0;
        % for each neighboring pixels of (i, j)
        for x = max(i-r, 1) : min(i+r, M)
            for y = max(j-r, 1): min(j+r, N)
                if (i==x)&&(j==y)
                    continue;
                end
                % if the neighbor is also unlabeled, then skip it
                if Omega(x, y) == 0
                    continue;
                end
                difference = abs(B(i, j) - B(x, y)) / (abs(B(x, y)) + 1e-2);
                distance = sqrt((x-i)^2 + (y-j)^2);
                if (difference < tau) 
                    isFindSimilarNeighbor = 1;
                    weight = exp(- difference^2 * distance);
                    G(i, j) = weight * D(x, y);
                    G(i, j+N) = weight * D(x, y+N);
                    G(i, j+2*N) = weight * D(x, y+2*N);
                end
            end
        end
        % if no similar labeled neighbor
        if isFindSimilarNeighbor == 0
            continue;
        end
        weight = B(i, j) * 3 / (G(i, j) + G(i, j+N) + G(i, j+2*N));
        G(i, j) = weight * G(i, j);
        G(i, j+N) = weight * G(i, j+N);
        G(i, j+2*N) = weight * G(i, j+2*N);
    end
end



end