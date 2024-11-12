%
clear
clc
close all
addpath('image/');
%
gamma  = 0.04; mu = 1/(8*gamma);
tol=1e-5; maxiter.out=100; maxiter.in = 10;

x_t  = double(imread('Monarch256.tif'))./255;
%x_t  = rgb2gray(x_t);
%name = {'lena'};
figure;imshow(x_t);title('‘≠Õº');
                   
[m,n]=size(x_t);
rng('default');
v1 = randn(m,n); v2 = randn(m,n);
v3 = randn(m,n); v4 = randn(m,n);
            
x_n = x_t+gamma*v3./v4;
x_t = min(max(0,x_t),1);
x_n = min(max(-10,x_n),20);
alpha = 0.15; M = 0.12; beta = 15; r = 45; 
figure;imshow(x_n);title('‘Î…˘Õº'); 

tic;
[u,g1,g2] = Cauchy_denoising(x_n,alpha,M,r,beta,tol,maxiter,sqrt(gamma),mu);  
time = toc;
u  = min(max(0,u),1);
ps = psnr_fun(u,x_t);
ss = ssim((u*255),(x_t*255));
figure;imshow(u);title('∏¥‘≠Õº');

display(sprintf('psnr=%.4f,ssim=%.4f,alpha=%.2f,M=%.2f,beta=%.0f,r=%.0f', ps, ss, alpha, M, beta, r))
display(sprintf('=================================='))
        

