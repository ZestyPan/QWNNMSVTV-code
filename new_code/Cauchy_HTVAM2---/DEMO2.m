%
clear
clc
close all
addpath('image/');
%
kernel = fspecial('gaussian', [9,9], 1);
gamma  = 0.02; mu = 1/(8*gamma);
tol=1e-5; maxiter.out=60; maxiter.in = 10;


x_t=double(imread('starfish256.tif'))./255;figure;imshow(x_t);title("‘≠Õº");
name = {'house'};
            
[m,n]=size(x_t);
rng('default');
v1 = randn(m,n); v2 = randn(m,n);
v3 = randn(m,n); v4 = randn(m,n);
            
x_n = imfilter(x_t,kernel,'circular','conv')+gamma*v3./v4;
figure;imshow(x_n);title("‘Î…˘Õº");

x_t = min(max(0,x_t),1);
x_n = min(max(-10,x_n),20);
M = 0.07; r = 5; alpha = 0.14; beta = 18;  

tic;
[u,g1,g2] = Cauchy_deblurring(x_n,kernel,alpha,M,r,beta,tol,maxiter,sqrt(gamma),mu);  
time = toc;
u  = min(max(0,u),1);figure;imshow(u);title("∏¥‘≠Õº");
ps = psnr_fun(u,x_t);
ss = ssim((u*255),(x_t*255));
            

  
 display(sprintf('psnr=%.4f,ssim=%.4f,alpha=%.2f,M=%.2f,beta=%.0f,r=%.0f', ps, ss, alpha, M, beta, r))
 display(sprintf('=================================='))
        

                              

