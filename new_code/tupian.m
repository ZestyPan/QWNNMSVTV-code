%% 读取待放大图像
% clc;clear;close all;
im = imread('QWNNM_SV_TV.jpg');
% im = imread('t2star_layer11_w_dy.bmp');

%% 在图像中局部显示矩形框
X    = 120; %矩形框左上角的横坐标   
Y    = 150; %矩形框左上角的纵坐标
dX   = 60;    
dY   = 60;  

bbox = [X, Y, dX, dY];
im_1 = insertShape(im, 'Rectangle', bbox, 'LineWidth', 3, 'Color', 'red');

%% 裁剪和局部并插值放大
scale     = 3; %1.8
im_crop   = imcrop(im, bbox); 
im_crop_b = imresize(im_crop, scale, 'bicubic'); 

%% 局部显示
[row_1, col_1, ~]                     = size(im);
[row_2, col_2, ~]                     = size(im_crop_b);
im_2                                  = im_1;

% 局部放大图放在左下角
%im_2(row_1-row_2+1:col_1, 1:col_2, :) = im_crop_b;
%bbox  = [1,row_1-row_2+1,col_2,row_2];

% 局部放大图放在右下角
im_2(row_1 - row_2 + 1:row_1, col_1 - col_2 + 1:col_1, :) = im_crop_b;
%im_2(row_1-row_2+1:col_1, col_1-col_2+1:row_1, :) = im_crop_b;
bbox  = [col_1-col_2+1,row_1-row_2+1,col_2,row_2];

%%
im_2  = insertShape(im_2, 'Rectangle', bbox, 'LineWidth', 4, 'Color', 'red');
figure;imshow(im_2);
% figure;imshow(im_crop_b);

%% 保存为.PNG格式文件
% imwrite(im_2,'layer14s_local.png');

