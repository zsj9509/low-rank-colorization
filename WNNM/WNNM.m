
function  [X] =  WNNM( Y, C, NSig, m, Iter )
    % [U,SigmaY,V] =   svd(full(Y),'econ');    
    [U,SigmaY,V] =   stablesvd(full(Y));
    PatNum       = size(Y,2);
    Temp         =   sqrt(max( diag(SigmaY).^2 - PatNum*NSig^2, 0 ));
    for i=1:Iter
        % Weight vector
        W_Vec    =   (C*sqrt(PatNum)*NSig^2)./( Temp + eps );               
        SigmaX   =  soft(SigmaY, diag(W_Vec));
        Temp     = diag(SigmaX);
    end
    X =  U*SigmaX*V' + m;     
return;
