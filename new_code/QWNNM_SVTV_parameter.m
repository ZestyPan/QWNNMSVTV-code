function[out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par)

%参数设置

epsi = 1*1e-6;%for stopping criteria 停止准则的阈值
deblur = 1;
[m,n,c] = size(imgO);
N = m*n*c;
r1     =       double2q(zeros(size(imgN))); %创建一个与噪声图像大小相同的全零矩阵，并将其转换为q格式
r1 = double2q(r1, 'inverse');


%% q子问题

    imgND = qQCTV_ADMM2(imgN, eta, beta1,beta2, iters, mu, deblur, H, epsi,r1,lamada);
    q = imgND;
    %ssim = ssim(imgND./255,imgO./255)
    %psnr = psnr(imgND./255,imgO./255)


%% g子问题
    q2double = @(X) double2q(X, 'inverse');%将量化的图像数据转换回双精度浮点数格式
    [m, n, c]  = size(imgO);%获取原始图像的尺寸
    TotalPatNum = (m-Par.patsize+1)*(n-Par.patsize+1);% 计算块的总数
    [Neighbor_arr, Num_arr, Self_arr] =	QWNNM_NeighborIndex(imgO, Par);%获取块的邻居索引
    NL_mat = zeros(Par.patnum,length(Num_arr));%创建一个与Num_arr长度相同的全零矩阵
    q       =       double2q(q);
    r1     =       double2q(zeros(size(imgN))); %创建一个与噪声图像大小相同的全零矩阵，并将其转换为q格式
    N_Img = q2double(q + r1./beta1);
    E_Img = N_Img;
    %E_Img = E_Img + Par.delta*(N_Img - E_Img);
    E_Img = E_Img + beta1*(N_Img - E_Img);
    [CurPat, Sigma_arr] = QWNNM_Im2Patch(E_Img, N_Img, Par);%将图像转换为块，并计算一些参数
    NL_mat = QWNNM_Block_matching(CurPat, Par, Neighbor_arr, Num_arr, Self_arr);%计算块之间的相似性
   
    if (mod(iters-1,Par.Innerloop)==0)%如果达到内部迭代次数，则执行以下操作
        Par.patnum = Par.patnum-10;                                             % Lower Noise level, less NL patches
        NL_mat = QWNNM_Block_matching(CurPat, Par, Neighbor_arr, Num_arr, Self_arr);%计算块之间的相似性

        % Caculate Non-local similar patches for each 
        %降低噪声水平，减少非局部相似补丁的数量，并计算每个补丁的非局部相似补丁
        if(iters==1)
            Sigma_arr = Par.nSig * ones(size(Sigma_arr));                       % First Iteration use the input noise Parameter
        end%如果是第一次迭代，则使用输入的噪声参数估计噪声方差。
    end      

    [EPat, W] = QWNNM_PatEstimation(NL_mat, Self_arr, Sigma_arr, CurPat, Par);%估计块的像素值
    E_Img = QWNNM_Patch2Im(EPat, W, Par.patsize, m, n, c);%将估计的块重构为图像
    %g = double2q(E_Img);  %g子问题
    g = E_Img;  

%% 拉格朗日乘子更新
        r1 = r1 + 0.001*(q - g);%拉格朗日乘子eta更新子问题

out = E_Img;
end

%% 用到的函数

% function imgND = qQCTV_ADMM2(imgN, eta, beta1,beta2, iters,mu, deblur, H, epsi,r1)
% 
% F = qctv2ctv(imgN,0);
% 
% [rows, cols, c] = size(F); 
% 
%     wx0 = zeros(rows, cols, c);
% 
%     wy0 = zeros(rows, cols, c);
%     
%     
%     lambdawx = zeros(rows, cols,c); 
%     lambdawy = zeros(rows, cols,c);
%     
%    wx = wx0;
%    wy = wy0;
%    
%    X_old = zeros(rows, cols,c); 
%    
% for iter = 1:iters
%     
%       X = qFFTsolutionofL2norm_QCTV(lambdawx,lambdawy,wx,wy,F,beta1,beta2,deblur,H,r1);
%       
% [DX,DY] = ForwardDM(X);
% 
% [wx,wy] = QCTV_shrinkage(DX,DY,lambdawx,lambdawy,eta,beta2,mu);
% 
% lambdawx = lambdawx+beta2*(wx-DX);
% lambdawy = lambdawy+beta2*(wy-DY);
% 
% 
%     if mod(iter, 2) == 0
%     
%         imgND = ctv2qctv(X,0);
%         imagesc(uint8(imgND*255));colormap gray;
%         drawnow
%     end
% 
% crit = norm(X(:,:,1)-X_old(:,:,1),2)/norm(X_old(:,:,1),2);
% 
%     if crit<epsi
%         
%          iter
%         break
%     end
% 
%     X_old = X;
% 
% end
% imgND = ctv2qctv(X,0);
% 
% end
% 
% 
% %% 
% function [vx,vy] = ForwardDM(v)
% 
% 
% vx = ForwardX(v);
% 
% vy = ForwardY(v);
% end
% 
% 
% function [dx] = ForwardX(v)
% 
% [Ny,Nx,Nc] = size(v);
% dx = zeros(size(v));
% dx(1:Ny-1,1:Nx-1,:)=( v(1:Ny-1,2:Nx,:) - v(1:Ny-1,1:Nx-1,:) );
% end
% 
% 
% 
% function [dy] = ForwardY(v)
% 
% [Ny,Nx,Nc] = size(v);
% dy = zeros(size(v));
% dy(1:Ny-1,1:Nx-1,:)=( v(2:Ny,1:Nx-1,:) - v(1:Ny-1,1:Nx-1,:) );
% end

