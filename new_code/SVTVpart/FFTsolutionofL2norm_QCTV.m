
function f = FFTsolutionofL2norm_QCTV(lambdawx,lambdawy,wx,wy,F,beta,deblur,h)


if deblur==0

for ch = 1:3
sizev = size(F(:,:,1));
otfDx = psf2otf([1,-1],sizev);
    otfDy = psf2otf([1;-1],sizev);
     
     Nomin = beta*(conj(otfDx).*fft2(wx(:,:,ch)+lambdawx(:,:,ch)/beta)+conj(otfDy).*fft2(wy(:,:,ch)+lambdawy(:,:,ch)/beta))+fft2(F(:,:,ch));
   
     Denom = beta*(abs(otfDx).^2 + abs(otfDy ).^2)+1;
   
     s = Nomin./Denom;
     f(:,:,ch) = real(ifft2(s));
end

else
    
 for ch = 1:3
sizev = size(F(:,:,1));
otfH = psf2otf(h,sizev);
otfDx = psf2otf([1,-1],sizev);
    otfDy = psf2otf([1;-1],sizev);
     
     Nomin = beta*(conj(otfDx).*fft2(wx(:,:,ch)+lambdawx(:,:,ch)/beta)+conj(otfDy).*fft2(wy(:,:,ch)+lambdawy(:,:,ch)/beta))+conj(otfH).*fft2(F(:,:,ch));
   
     Denom = beta*(abs(otfDx).^2 + abs(otfDy).^2)+abs(otfH).^2;
   
     s = Nomin./Denom;
     f(:,:,ch) = real(ifft2(s));
 end   
  
 
end









