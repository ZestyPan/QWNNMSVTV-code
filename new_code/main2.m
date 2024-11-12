%% 引入图像并加噪加模糊
% clc
% clear all
% close all
warning('off')

%for i=10:27
%load(strcat('D:\Zesty\quaternion\QWNNM_SV_TV\data\Set27_35_GB(15,1)\Set27_35_GB(15,1)',num2str(i),'.mat'))
%load(strcat('D:\Zesty\quaternion\QWNNM_SV_TV\data\Set27_35_MB(5,15)\Set27_35_MB(5,15)',num2str(i),'.mat'))
%load(strcat('D:\Zesty\quaternion\QWNNM_SV_TV\data\Set27_35_AB(3)\Set27_35_AB(3)',num2str(i),'.mat'))

load('Set27_(15,1)2.mat');
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

Par.nSig      =   nSig; 
Par.c         =   1.3*sqrt(2); 
Par.Innerloop           =   1;%迭代次数
Par.ReWeiIter           =   1;
Par.Iter                =   5;%迭代次数
Par.SearchWin           =   30;                                   %非局部块搜索窗口的大小 Non-local patch searching window
Par.patsize             =    7;%块的尺寸
Par.patnum              =       330; %每个块的数量     
Par.lamada        =       0.28;%110
lamada        =       0.28;%110
Par.step      =  2;

