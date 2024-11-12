function X = qFFTsolutionofL2norm_QCTV(lambdawx,lambdawy,wx,wy,F,beta1,beta2,deblur,H,r1,lamada)
% g = F
if deblur==0

for ch = 1:3
sizev = size(F(:,:,1));%获取输入图像的大小
    otfDx = psf2otf([1,-1],sizev);
    otfDy = psf2otf([1;-1],sizev);%将一维的差分算子转换为频率域的算子，用于计算梯度

    %% deblur==0 最原始的代码
%     Nomin = beta*(conj(otfDx).*fft2(wx(:,:,ch)+lambdawx(:,:,ch)/beta)+conj(otfDy).*fft2(wy(:,:,ch)+lambdawy(:,:,ch)/beta))+fft2(F(:,:,ch));
%     Denom = beta*(abs(otfDx).^2 + abs(otfDy ).^2)+1;

%能运行处理且结果较好的代码
    Nomin =beta1*fft2(F(:,:,ch))-r1(:,:,ch)+beta2*(conj(otfDx).*fft2(wx(:,:,ch)+lambdawx(:,:,ch)/beta2)+conj(otfDy).*fft2(wy(:,:,ch)+lambdawy(:,:,ch)/beta2))+lamada.*fft2(F(:,:,ch));
    Denom = beta1+beta2*(abs(otfDx).^2 + abs(otfDy ).^2);
    %Denom = lamada+beta1+beta2*(abs(otfDx).^2 + abs(otfDy ).^2);
     s = Nomin./Denom;
     f(:,:,ch) = real(ifft2(s));
end

else
    
 for ch = 1:3
sizev = size(F(:,:,1));
otfH = psf2otf(H,sizev);
otfDx = psf2otf([1,-1],sizev);
otfDy = psf2otf([1;-1],sizev);
     

%% deblur==1
%最原始的代码
%     Nomin = beta*(conj(otfDx).*fft2(wx(:,:,ch)+lambdawx(:,:,ch)/beta)+conj(otfDy).*fft2(wy(:,:,ch)+lambdawy(:,:,ch)/beta))+conj(otfH).*fft2(F(:,:,ch));
%     Denom = beta*(abs(otfDx).^2 + abs(otfDy).^2)+abs(otfH).^2;
   
% 能运行 & 结果较好
%     Nomin = beta2*(conj(otfDx).*fft2(wx(:,:,ch)+lambdawx(:,:,ch)/beta2)+conj(otfDy).*fft2(wy(:,:,ch)+lambdawy(:,:,ch)/beta2))+conj(otfH).*fft2(F(:,:,ch));
%     Denom = beta2*(abs(otfDx).^2 + abs(otfDy).^2)+abs(otfH).^2 ;

% q子问题的公式
    Nomin =beta1*fft2(F(:,:,ch))-r1(:,:,ch)+beta2*(conj(otfDx).*fft2(wx(:,:,ch)+lambdawx(:,:,ch)/beta2)+conj(otfDy).*fft2(wy(:,:,ch)+lambdawy(:,:,ch)/beta2))+lamada.*conj(otfH).*fft2(F(:,:,ch));
    %Nomin       =       real(Nomin);
    Denom = beta1+beta2*(abs(otfDx).^2 + abs(otfDy).^2)+lamada.*abs(otfH).^2 ;
     s = Nomin./Denom;
     f(:,:,ch) = real(ifft2(s));
     %f=real(ifft2(s));
 end   

end
X=f;
end