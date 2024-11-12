function  [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada)
% Par.nSig      =   nSig;                                 % Variance of the noise image
% Par.SearchWin =   30;                                   % Non-local patch searching window
% Par.c       =   1.3*sqrt(2);%%sqrt(2)                            % Constant num for the weight vector
% Par.Innerloop =   2;      % 块匹配之间的内部循环次数InnerLoop Num of between re-blockmatching
% Par.ReWeiIter =   1;%%3   %加权迭代次数



%Par.c       =   c;%%sqrt(2)                            % Constant num for the weight vector
Par.c         =   1.3*sqrt(2);%%sqrt(2)
Par.Innerloop =   1;      % 块匹配之间的内部循环次数InnerLoop Num of between re-blockmatching
Par.ReWeiIter =   1;%%3   %加权迭代次数
Par.Iter      =  5;


Par.nSig      =   nSig;                                 % Variance of the noise image
Par.SearchWin =   SearchWin;
Par.patsize       =   patsize;
Par.patnum        =   patnum; 
Par.lamada        =   lamada;
Par.step         = step;


% if nSig<=20
%     %Par.patsize       =   6;                            % Patch size
%     %Par.patnum        =   80;     %70                      % Initial Non-local Patch number
%     %Par.Iter          =   9;          %8                  % total iter numbers
%     %Par.lamada        =   0.54;                         % Noise estimete Parameter
%     Par.patsize       =   patsize;
%     Par.patnum        =   patnum;
%     Par.Iter          =   9;    
%     Par.lamada        =   lamada;
%     Par.step         = step;
% 
% elseif nSig <= 40
% %     Par.patsize       =   7;
% %     Par.patnum        =   90;
% %     Par.Iter          =   12;
% %     Par.lamada        =   0.56; 
%     Par.patsize       =   patsize;
%     Par.patnum        =   patnum;
%     Par.Iter          =   12;    
%     Par.lamada        =   lamada;
%     Par.step         = step;
% elseif nSig<=60
% %     Par.patsize       =   8;
% %     Par.patnum        =   120;
% %     Par.Iter          =   14;
% %     Par.lamada        =   0.58; 
%     Par.patsize       =   patsize;
%     Par.patnum        =   patnum;
%     Par.Iter          =   14;    
%     Par.lamada        =   lamada;
%     Par.step         = step;
% else
% %     Par.patsize       =   9;
% %     Par.patnum        =   140;
% %     Par.Iter          =   14;
% %     Par.lamada        =   0.58; 
%     Par.patsize       =   patsize;
%     Par.patnum        =   patnum;
%     Par.Iter          =   14;    
%     Par.lamada        =   lamada;
%     Par.step         = step;
% end

% Par.step      =   floor((Par.patsize)/2-1);                   
