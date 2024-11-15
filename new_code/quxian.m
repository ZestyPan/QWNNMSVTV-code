% warning('off')
% 
% load('Set27_(15,1)2.mat');
% Sig   = 25;% 15
% H     = fspecial('gaussian',15,1); 
% 
% %% 参数+初始化
% [m,n,c] = size(imgO);
% N = m*n*c;
% tau = 0.0010;
% eta = tau*sqrt((Sig^2)*N);
% tau=[0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.0010]
% x = [2.1651, 4.3301, 6.4952, 8.6603, 10.8253, 12.9904, 15.1554, 17.3205,  19.4856, 21.6506]; % x轴上的点
% y = [27.6008, 27.6188, 27.6288, 27.6293, 27.6291, 27.6289, 27.6268, 27.6270, 27.6270, 27.6270]; % y轴上的点
% 
% plot(x,y,'-+r');
% 
% 
% xlabel(''); % 设置x轴标签
% ylabel('PSNR'); % 设置y轴标签


%x = [0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40]; %lambda
%x = [2,3,4,5,6,7,8,9,10]  %patsize
%x = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]; %eta
x = [120,150,180,210,240,270,300,330,360,390,420,450,480,510,540];

%img1
y = [23.0799,23.3924,23.6032,23.9147,24.1019,24.1272,24.1307,24.1321,24.1315,24.1290,24.1245,24.0183,23.8153,23.4052,23.0957];
plot(x,y,'r','linewidth',2);
hold on;
%img2
y = [26.6451,26.9726,27.2923,27.4203,27.5273,27.6311,27.6329,27.6332,27.6289,27.6221,27.6162,27.5978,27.5681,27.5252,27.4790];
plot(x,y,'m','linewidth',2);
hold on;
%img4
y = [28.0982,28.3364,28.6629,28.9782,29.0896,29.1962,29.2051,29.2132,29.2228,29.2285,29.2272,29.2257,29.1291,29.1008,28.9376];
plot(x,y,'g','linewidth',2);
hold on;
%img5
y = [32.2120 ,32.5563,32.7822,32.9892,33.0948,33.1131,33.1278,33.1376,33.1345,33.1254,33.1159 ,33.0921 ,32.9041 ,32.7494,32.6602];
plot(x,y,'b','linewidth',2);
hold on;
%img6
y = [25.0698,25.2965,25.5302,25.7592 ,25.8894 ,25.8951,25.9238,25.9224,25.9141,25.9002,25.8951,25.7657,25.6365,25.4998,25.2678 ];
plot(x,y,'k','linewidth',2);
hold on;
%img13
y = [25.8978,26.0848,26.1169,26.2346,26.3480 ,26.3575,26.3629,26.3645,26.3733,26.3684,26.3619,26.2535,26.1421,26.0270 ,25.9957];
plot(x,y,'c','linewidth',2);
hold on;







xlabel(''); % 设置x轴标签
ylabel('PSNR'); % 设置y轴标签
legend('img1','img2','img4','img5','img6','img13');