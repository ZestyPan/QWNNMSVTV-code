%
clear
clc
close all
% addpath('image/');
%
gamma  = 5;
% mu = 1/(8*gamma);
mu = 10;
tol=1e-5; maxiter.out=100; maxiter.in = 10;

% x_t  = double(imread('starfish256.tif'))./255;
% name = {'cameraman'};
%                   
% [m,n]=size(x_t);
% rng('default');
% v1 = randn(m,n); v2 = randn(m,n);
% v3 = randn(m,n); v4 = randn(m,n);
%             
% x_n = x_t+gamma*v3./v4;
% x_t = min(max(0,x_t),1);
% x_n = min(max(-10,x_n),20);

%load('baby_01','O_Img','Y');
load('Set27_01.mat');

x_n = max(0,min(imgN,255));
alpha = 0.14; M = 0.12; beta = 16; r = 39; 

% for ImgNo = 2
%     switch ImgNo
        %case 1
lambda =1;h = 1;max =10;data=[];
while lambda < max
    fprintf( 'Estimated Image: lambda = %3f\n',lambda);
    t1 = cputime;
     %[out] = OGSTV_denoising( Y, r, lambda, mu, beta1, beta2, K ); 
     [out] = QWNNM_SVTV_parameter( Y, r, lambda, mu, beta1, beta2, K ); 
    t1 = cputime - t1; 
    s1 = psnr(uint8(out.sol),uint8(O_Img));
    fprintf('pnsr = %.2f\n',s1);
    data=[data;lambda,s1];
    lambda = lambda +h;
end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');
%         case 2
% beta1 =500;h = 100;max =1500;data=[];
% while beta1 < max
%     fprintf( 'Estimated Image: beta1 = %3f\n',beta1);
%     t1 = cputime;
%  [out] = OGSTV_denoising( Y, r, lambda, mu, beta1, beta2, K );
%     t1 = cputime - t1; 
%     s1 =psnr(uint8(out.sol),uint8(O_Img));
%     fprintf('pnsr = %.2f\n',s1);
%     data=[data;beta1,s1];
%     beta1 = beta1 +h;
% end
% figure(1);
% plot(data(:,1),data(:,2),'-+r');
% ylabel('PSNR');
%         case 3
% beta2 =1000;h = 100;max =2500;data=[];
% while beta2 < max
%     fprintf( 'Estimated Image: beta2 = %3f\n',beta2);
%     t1 = cputime;
%    [out] = OGSTV_denoising( Y, r, lambda, mu, beta1, beta2, K );
%     t1 = cputime - t1; 
%     s1 = psnr(uint8(out.sol),uint8(O_Img));
%     fprintf('pnsr = %.2f\n',s1);
%     data=[data;beta2,s1];
%     beta2 = beta2 +h;
% end
% figure(1);
% plot(data(:,1),data(:,2),'-+r');
% ylabel('PSNR');
%         case 4
% mu =5;h = 1;max =25;data=[];
% while mu < max
%     fprintf( 'Estimated Image: mu = %3f\n',mu);
%     t1 = cputime;
%    [out] = OGSTV_denoising( Y, r, lambda, mu, beta1, beta2, K );
%     t1 = cputime - t1; 
%     s1 = psnr(uint8(out.sol),uint8(O_Img));
%     fprintf('pnsr = %.2f\n',s1);
%     data=[data;mu,s1];
%     mu = mu +h;
% end
% figure(1);
% plot(data(:,1),data(:,2),'-+r');
% ylabel('PSNR');
%     end
% end

% tic;
% [u,g1,g2] = Cauchy_denoising(x_n,alpha,M,r,beta,tol,maxiter,gamma,mu);  
% time = toc;

% u  = min(max(0,u),1);
ps = psnr(uint8(u),uint8(O_Img));
ss = ssim(uint8(u),uint8(O_Img));
% u1=u(90:153,120:183);
% u1=imresize(u1,4);
% imwrite(u1,'straw-2yang-big.png','png');
% imwrite(u,'straw-2yang.png','png');


display(sprintf('psnr=%.4f,ssim=%.4f,alpha=%.2f,M=%.2f,beta=%.0f,r=%.0f', ps, ss, alpha, M, beta, r))
display(sprintf('=================================='))
        

