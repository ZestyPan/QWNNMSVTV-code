function IND = QCTV_ADMM2(IN, alpha, beta, itera,mu, deblur, h, epsi)

F = qctv2ctv(IN,0);

[rows, cols, c] = size(F); 

    wx0 = zeros(rows, cols, c);

    wy0 = zeros(rows, cols, c);
    
    
    lambdawx = zeros(rows, cols,c); 
    lambdawy = zeros(rows, cols,c);
    
   wx = wx0;
   wy = wy0;
   
   X_old = zeros(rows, cols,c); 
   
for iter = 1:itera
    
      X = FFTsolutionofL2norm_QCTV(lambdawx,lambdawy,wx,wy,F,beta,deblur,h);
      
[DX,DY] = ForwardDM(X);

[wx,wy] = QCTV_shrinkage(DX,DY,lambdawx,lambdawy,alpha,beta,mu);

lambdawx = lambdawx+beta*(wx-DX);
lambdawy = lambdawy+beta*(wy-DY);


if mod(iter, 2) == 0

    IND = ctv2qctv(X,0);
    imagesc(uint8(IND*255));colormap gray;
    drawnow
end

crit = norm(X(:,:,1)-X_old(:,:,1),2)/norm(X_old(:,:,1),2);

if crit<epsi
    
iter
    break
end

X_old = X;

end
IND = ctv2qctv(X,0);


function [vx,vy] = ForwardDM(v)


vx = ForwardX(v);

vy = ForwardY(v);


function [dx] = ForwardX(v)

[Ny,Nx,Nc] = size(v);
dx = zeros(size(v));
dx(1:Ny-1,1:Nx-1,:)=( v(1:Ny-1,2:Nx,:) - v(1:Ny-1,1:Nx-1,:) );



function [dy] = ForwardY(v)

[Ny,Nx,Nc] = size(v);
dy = zeros(size(v));
dy(1:Ny-1,1:Nx-1,:)=( v(2:Ny,1:Nx-1,:) - v(1:Ny-1,1:Nx-1,:) );




function [wx,wy] = QCTV_shrinkage(DX,DY,lambdawx,lambdawy,alpha,beta,mu)


s1 = sqrt((DX(:,:,1)-lambdawx(:,:,1)/beta).^2+(DX(:,:,2)-lambdawx(:,:,2)/beta).^2+(DY(:,:,1)-lambdawy(:,:,1)/beta).^2+(DY(:,:,2)-lambdawy(:,:,2)/beta).^2);
s1(s1 == 0) = 1e-3;

s2 = sqrt((DX(:,:,3)-lambdawx(:,:,3)/beta).^2+(DY(:,:,3)-lambdawy(:,:,3)/beta).^2);
s2(s2 == 0) = 1e-3;

wx(:,:,3) = max(0, s2-alpha*mu/beta).*(DX(:,:,3)-lambdawx(:,:,3)/beta)./s2;
wx(:,:,1) = max(0, s1-alpha/beta).*(DX(:,:,1)-lambdawx(:,:,1)/beta)./s1;
wx(:,:,2) = max(0, s1-alpha/beta).*(DX(:,:,2)-lambdawx(:,:,2)/beta)./s1;

wy(:,:,3) = max(0, s2-alpha*mu/beta).*(DY(:,:,3)-lambdawy(:,:,3)/beta)./s2;
wy(:,:,1) = max(0, s1-alpha/beta).*(DY(:,:,1)-lambdawy(:,:,1)/beta)./s1;
wy(:,:,2) = max(0, s1-alpha/beta).*(DY(:,:,2)-lambdawy(:,:,2)/beta)./s1;














