function  [par]=QWNNM_ParSet(nSig)

par.nSig      =   nSig;                                 % Variance of the noise image
par.SearchWin =   30;                                   %非局部块搜索窗口的大小 Non-local patch searching window
par.delta     =   0.1;    %0.1  设置参数 delta，在每次迭代中用于调整阈值 % Parameter between each iter
par.c         =   1.3*sqrt(2);%%sqrt(2)  设置常数 c，用于计算权重向量% Constant num for the weight vector
par.Innerloop =   1;    %设置内部迭代次数，即重新进行块匹配的次数 InnerLoop Num of between re-blockmatching
par.ReWeiIter =   1;%%3     重新计算权重的迭代次数
if nSig<=20  %根据不同的噪声水平设置不同的参数值
    par.patsize       =   6;  % 噪声水平设置块的大小Patch size
    par.patnum        =   80;     %70  设置初始的非局部块数量 % Initial Non-local Patch number
    par.Iter          =   9;          %8 设置总的迭代次数% total iter numbers
    par.lamada        =   0.54;                         % Noise estimete parameter
elseif nSig <= 40
    par.patsize       =   7;
    par.patnum        =   90;
    par.Iter          =   12;
    par.lamada        =   0.56; 
elseif nSig<=60
    par.patsize       =   8;
    par.patnum        =   120;
    par.Iter          =   14;
    par.lamada        =   0.58; 
else
    par.patsize       =   9;
    par.patnum        =   140;
    par.Iter          =   14;
    par.lamada        =   0.58; 
end

par.step      =   floor((par.patsize)/2-1);       %计算块匹配时的步长            
