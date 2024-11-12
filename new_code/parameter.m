%
% clear
% clc
% close all
warning('off')

load('Set27_(15,1)6.mat');
imgO  = double(imread('6.png'));
H     = fspecial('gaussian',15,1);  %创建一个高斯模糊

[m,n,c] = size(imgO);
N = m*n*c;
tol=1e-5; maxiter.out=100; maxiter.in = 10;
Sig=25;
theta = 8.5;
nSig = Sig*theta;
mu = 0.1; 
beta1 = 7.5; %QWNNM过程中
beta2 = 0.0004; %SV-TV过程中
iters = 100;%iteration numbers迭代次数
Iter        =   5;%迭代次数
%主要要调的系数参数
lamada              =       0.28;%正则化参数      
tau = 0.0004;
eta = tau*sqrt((Sig^2)*N);

%要调的块的参数
SearchWin   =   30;                                   % Non-local patch searching window
patsize     =   7;                            % Patch size
patnum      =   330;     %70                      % Initial Non-local Patch number
step        =   2;


for ImgNo = 8
     switch ImgNo
% %1：lamada 2：beta1 3：beta2 4：mu %5：eta 6：SearchWin 7：patsize 8：patnum 9：step
% 主要调1、5、6、7、8、9
case 1

lamada =0.21;h = 0.01;max =0.33;data=[];
    while lamada < max
        fprintf( 'Estimated Image: lamada = %3f\n',lamada);
        tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;lamada,s1];
        lamada = lamada +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');

case 2
beta1 =0;h = 1;max =10;data=[];
    while beta1 < max
        fprintf( 'Estimated Image: beta1 = %3f\n',beta1);
        tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;beta1,s1];
        beta1 = beta1 +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');

case 3
beta2 =1000;h = 100;max =2500;data=[];
    while beta2 < max
        fprintf( 'Estimated Image: beta2 = %3f\n',beta2);
        tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;beta2,s1];
        beta2 = beta2 +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');

case 4
mu =5;h = 1;max =25;data=[];
    while mu < max
        fprintf( 'Estimated Image: mu = %3f\n',mu);
        tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;mu,s1];
        mu = mu +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');

case 5
%tau =0.0001;h = 0.0001;max =0.0010;data=[];
%eta = tau*sqrt((Sig^2)*N);
%eta = 2.1651;h = 2.1651;max = 22;data=[];
eta = 3.6687;h = 3.6687;max = 37;data=[];
    while eta < max
        fprintf( 'Estimated Image: eta = %3f\n',eta);
       tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;eta,s1];
        eta = eta +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');

%%% 以下是块的参数设置

case 6
SearchWin =0;h = 1;max =10;data=[];
    while SearchWin < max
        fprintf( 'Estimated Image: SearchWin = %3f\n',SearchWin);

       tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;SearchWin,s1];
        SearchWin = SearchWin +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');


case 7
patsize =2;h = 1;max =12;data=[];
    while patsize < max
        fprintf( 'Estimated Image: patsize = %3f\n',patsize);

        tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;patsize,s1];
        patsize = patsize +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');

case 8
patnum =210;h = 30;max =460;data=[];
    while patnum < max
        fprintf( 'Estimated Image: patnum = %3f\n',patnum);
        tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;patnum,s1];
        patnum = patnum +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');


case 9
step =0;h = 1;max =10;data=[];
    while step < max
        fprintf( 'Estimated Image: step = %3f\n',step);
        tic
        t1 = cputime;
        [Par]=QWNNM_SVTV_ParSet(nSig,SearchWin,patsize,patnum,step,lamada);        
        [out] = QWNNM_SVTV_parameter(imgO,imgN,H,lamada, eta , mu, beta1, beta2,iters,Par);
        t1 = cputime - t1; 
        toc
        s1 = psnr(out./255,imgO./255);
        s2 = ssim(out./255,imgO./255);
        fprintf('pnsr = %.4f，ssim = %.4f\n',s1,s2);
        data=[data;step,s1];
        step = step +h;
    end
figure(1);
plot(data(:,1),data(:,2),'-+r');
ylabel('PSNR');
    end
end



