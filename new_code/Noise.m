
%  O_Img = double(imread('Einstein.tiff'));
% gamma     = 5;
% Generate Cauchy noisy image.
% rng('default');
% eta1 = randn(size(O_Img));
% eta2 = randn(size(O_Img));
% noise = gamma.*eta1./eta2;
% Y = O_Img + noise ;

Sig   = 35;
imgO  = double(imread('5.png'));
%H     = fspecial('gaussian',15,1);  %创建一个高斯模糊
%H = fspecial('motion',5,15); % 创建一个运动模糊
H = fspecial('average',3); % 创建一个平均模糊
[imgB, H_full] = addblur(H, imgO);%将原始图像进行模糊处理，得到模糊图像imgB和模糊核H_full
imgN  = imgB + Sig * randn(size(imgB));% 模糊图像再加上高斯噪声 
% save('Einstein_01.mat','O_Img','Y');
%save('Set27_35_AB(3)27.mat','imgO','imgN');
figure;imshow(uint8(imgO));
figure;imshow(uint8(imgN));