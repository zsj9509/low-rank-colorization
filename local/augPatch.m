function [AugPat] = augPatch(pat, pos, imSize, c)
    m = imSize(1);
    n = imSize(2);
    
    pos(1,:) = c*pos(1,:)/n;
    pos(2,:) = c*pos(2,:)/m;

    pat = single(pat);
    pos = single(pos);

    AugPat = cat(1, pat, pos);
end