function imgND = qQCTV_ADMM2(imgN, eta, beta1,beta2, iters,mu, deblur, H, epsi,r1,lamada)

F = qctv2ctv(imgN,0);

[rows, cols, c] = size(F); 

    wx0 = zeros(rows, cols, c);

    wy0 = zeros(rows, cols, c);
    
    
    lambdawx = zeros(rows, cols,c); 
    lambdawy = zeros(rows, cols,c);
    
   wx = wx0;
   wy = wy0;
   
   X_old = zeros(rows, cols,c); 
   
for iter = 1:iters
    
     % X = qFFTsolutionofL2norm_QCTV(lambdawx,lambdawy,wx,wy,F,beta1,beta2,deblur,H,r1);
      X = qFFTsolutionofL2norm_QCTV(lambdawx,lambdawy,wx,wy,F,beta1,beta2,deblur,H,r1,lamada);
      
[DX,DY] = ForwardDM(X);

[wx,wy] = QCTV_shrinkage(DX,DY,lambdawx,lambdawy,eta,beta2,mu);

lambdawx = lambdawx+beta2*(wx-DX);
lambdawy = lambdawy+beta2*(wy-DY);


    if mod(iter, 2) == 0
    
        imgND = ctv2qctv(X,0);
        imagesc(uint8(imgND*255));colormap gray;
        drawnow;
    end

crit = norm(X(:,:,1)-X_old(:,:,1),2)/norm(X_old(:,:,1),2);

    if crit<epsi
        
         iter;
        break
    end

    X_old = X;

end
imgND = ctv2qctv(X,0);

end


%% 
function [vx,vy] = ForwardDM(v)


vx = ForwardX(v);

vy = ForwardY(v);
end


function [dx] = ForwardX(v)

[Ny,Nx,Nc] = size(v);
dx = zeros(size(v));
dx(1:Ny-1,1:Nx-1,:)=( v(1:Ny-1,2:Nx,:) - v(1:Ny-1,1:Nx-1,:) );
end



function [dy] = ForwardY(v)

[Ny,Nx,Nc] = size(v);
dy = zeros(size(v));
dy(1:Ny-1,1:Nx-1,:)=( v(2:Ny,1:Nx-1,:) - v(1:Ny-1,1:Nx-1,:) );
end
