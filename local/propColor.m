function [ Omega, Obvs] = propColor(gImg, Omega, Obvs, tol)
% gImg:  [m, n] gray image
% Omega: [m, n]

[m, n] = size(gImg);
if(numel(Omega) > numel(gImg))
    Omega = Omega(:, 1:n);
end
Omega = double(Omega > 0);

fprintf('before %2.2f \n', nnz(Omega)/numel(Omega));

[row, col] = find(Omega == 1);
for i = 1:numel(row)
    c = col(i);
    r = row(i);
    ctr = gImg(r, c);
    obs = Obvs(r, c, :);
    
    cc = c - 1;
    cr = r - 1;
    if(cr > 0 && cc > 0)
        lt = gImg(cr, cc);
        err = (ctr - lt)^2;
        if(err < tol && Omega(cr, cc) ~= 1)
            Omega(cr, cc) = 1;
            Obvs(cr, cc, :) = obs;
        end
    end
    if(r + 1 < m + 1 && c - 1 > 0)
        lb = gImg(r + 1, c - 1);
        err = (ctr - lb)^2;
        if(err < tol && Omega(r + 1, c - 1) ~= 1)
            Omega(r + 1, c - 1) = 1;
            Obvs(r + 1, c - 1, :) = obs;
        end
    end
    if(r + 1 < m + 1 && c + 1 < n + 1 )
        rb = gImg(r + 1, c + 1);
        err = (ctr - rb)^2;
        if(err < tol && Omega(r + 1, c + 1) ~= 1)
            Omega(r + 1, c + 1) = 1;
            Obvs(r + 1, c + 1, :) = obs;
        end
    end
    if(r - 1 > 0 && c + 1 < n + 1)
        rt = gImg(r - 1, c + 1);
        err = (ctr - rt)^2;
        if(err < tol && Omega(r - 1, c + 1) ~= 1)
            Omega(r - 1, c + 1) = 1;
            Obvs(r - 1, c + 1, :) = obs;
        end
    end
end

fprintf('after %2.2f \n', nnz(Omega)/numel(Omega));

end


