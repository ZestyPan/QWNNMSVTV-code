%贾志刚 SV-TV

clc
clear all
close all


%============================================
%load data
I = double(imread('10.png'))./255;

load 10n05
%=============================================


%==============================================
%parameter setting
[m,n,c] = size(I);
N = m*n*c;

deblur = 0;
h = fspecial('gaussian',[15,15],1.);%blur kernal

sigma = 0.05;%noise level

mu = 0.1;%mu for the inside value correction term

tau = 0.01;
alpha = tau*sqrt((sigma^2)*N);%regularization parameter

beta = 0.01;%for ADMM penalty

iters = 100;%iteration numbers

epsi = 1*1e-6;%for stopping criteria

%===============================================

%===============================================
%SVTV restoration
IND = QCTV_ADMM2(ID, alpha, beta, iters, mu, deblur, h, epsi);

ssim = ssim(IND,I)
psnr = psnr(IND,I)

%===============================================










