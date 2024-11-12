%E_Temp = WNNM(Temp, Par.c, Sigma_arr(Self_arr(i)), M_Temp, Par.ReWeiIter);

% function  [X] =  WNNM( Y, C, NSig, m, Iter )
%     [U,SigmaY,V] =   svd(full(Y),'econ');    
%     PatNum       = size(Y,2);
%     TempC  = C*sqrt(PatNum)*2*NSig^2;
%     [SigmaX,svp] = ClosedWNNM(SigmaY,TempC,eps);                        
%     X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)' + m;     
% return;

function  [X] =  WNNM( Y, C, NSig, m, Iter ,E_Img)
    [U,SigmaY,V] =   svd(full(E_Img),'econ');    
    PatNum       = size(Y,2);
    TempC  = C*sqrt(PatNum)*2*NSig^2;
    [SigmaX,svp] = ClosedWNNM(SigmaY,TempC,eps);                        
    X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';   
return;

