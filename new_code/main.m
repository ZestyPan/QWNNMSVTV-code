%% 引入图像并加噪加模糊
% clc
% clear all
% close all
warning('off')

%for i=10:27
%load(strcat('D:\Zesty\quaternion\QWNNM_SV_TV\data\Set27_35_GB(15,1)\Set27_35_GB(15,1)',num2str(i),'.mat'))
%load(strcat('D:\Zesty\quaternion\QWNNM_SV_TV\data\Set27_35_MB(5,15)\Set27_35_MB(5,15)',num2str(i),'.mat'))
%load(strcat('D:\Zesty\quaternion\QWNNM_SV_TV\data\Set27_35_AB(3)\Set27_35_AB(3)',num2str(i),'.mat'))

load('Set27_(15,1)13.mat');
Sig   = 25;% 15
%H     = fspecial('gaussian',15,1);  %创建一个高斯模糊

% load('Set27_(25,1.6)1.mat' );
H     = fspecial('gaussian',15,1); 
%H = fspecial('motion',5,15); % 创建一个运动模糊
%H = fspecial('average',3);  % 创建一个平均模糊

 %imgO  = double(imread('2.png'));
% [imgB, H_full] = addblur(H, imgO);%将原始图像进行模糊处理，得到模糊图像imgB和模糊核H_full
% imgN  = imgB + Sig * randn(size(imgB));% 模糊图像再加上高斯噪声 
PSNR  = psnr(imgN./255, imgO./255); %计算观测图像imgB的峰值信噪比
SSIM  = ssim(imgN./255, imgO./255); %计算观测图像imgB的峰值信噪比

%fprintf( 'i = %d\n', i);
fprintf( 'Noisy Image: Sig = %2.4f, PSNR = %2.4f ,SSIM = %2.4f \n\n', Sig, PSNR ,SSIM);
figure,imshow(uint8(imgO));
figure,imshow(uint8(imgN));
%figure;imshow(cat(2,uint8(imgO),uint8(imgN)));

%%% 观测图像(噪声+模糊)：imgN 清晰图像：imgO  模糊核：H  噪声：Sig   去噪后的图像：imgND

%% 参数+初始化
[m,n,c] = size(imgO);
N = m*n*c;
tau = 0.0004;
eta = tau*sqrt((Sig^2)*N);
beta1 = 8.5; %QWNNM过程中
beta2 = 0.0006; %SV-TV过程中
deblur = 1;
mu = 0.1;
epsi = 1*1e-6;%for stopping criteria 停止准则的阈值
iters = 100;%iteration numbers迭代次数
theta = 8.5;
r1     =       double2q(zeros(size(imgN))); %创建一个与噪声图像大小相同的全零矩阵，并将其转换为q格式
r1 = double2q(r1, 'inverse');
%%% beta1：QWNNM中的惩罚参数；beta2：SV-TV中的惩罚参数；mu：wx，wy子问题中的阿尔法
%%% eta：SV-TV前的稀疏权重 相当于SV-TV代码中的alpha
%%% Par.lamada是最开始的系数lamada；r1是变量g的拉格朗日乘子

nSig = Sig*theta;   %计算噪声标准差 

% SearchWin = 30;
% patsize   =  7;
% patnum    =  180;
% step      =  2;
% lamada    =  153;%正则化参数
% Par   = QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);%置QWNNM算法相关的参数。    

%Par   = QWNNM_ParSet(nSig);%置QWNNM算法相关的参数。
Par.nSig      =   nSig; 
Par.c         =   1.3*sqrt(2); 
Par.Innerloop           =   1;%迭代次数
Par.ReWeiIter           =   1;
Par.Iter                =   5;%迭代次数
Par.SearchWin           =   30;                                   %非局部块搜索窗口的大小 Non-local patch searching window
Par.patsize             =    7;%块的尺寸
Par.patnum              =       480; %每个块的数量     
Par.lamada        =       0.28;%110
lamada        =       0.28;%110
Par.step      =  2;

tic
%for iter = 1:8%开始迭代
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
       % r1 = r1 + 0.001*(q - g);%拉格朗日乘子eta更新子问题

%end
toc
%% 输出
% E_Img_rgb = repmat(E_Img, [1, 1, 3]);
% E_Img_cropped = E_Img_rgb(:, 1:500, :);
PSNR = psnr(E_Img./255,imgO./255);
SSIM=ssim(E_Img./255,imgO./255);
figure;imshow(cat(2,uint8(imgO),uint8(imgN),uint8(E_Img)));
drawnow;
fprintf( 'PSNR = %2.4f ,SSIM = %2.4f \n\n', PSNR ,SSIM);
fprintf( 'lamada = %2.4f,tau = %.4f\n', lamada, tau);
fprintf( 'Par.SearchWin = %2.4f ,Par.patsize = %2.4f,Par.patnum = %2.4f,Par.step = %.4f\n', Par.SearchWin ,Par.patsize, Par.patnum, Par.step);
title(['Original','               ','degraded             ',num2str(PSNR,'%2.2f'),'dB'],'FontSize',16);
figure,imshow(uint8(E_Img));
%end
