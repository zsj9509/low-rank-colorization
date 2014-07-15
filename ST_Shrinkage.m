function [Z] = ST_Shrinkage(Y, W)

Z = sign(Y).* max(0, abs(Y)-W) ;

end