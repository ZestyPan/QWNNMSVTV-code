function [wx,wy] = QCTV_shrinkage(DX,DY,lambdawx,lambdawy,eta,beta2,mu)


s1 = sqrt((DX(:,:,1)-lambdawx(:,:,1)/beta2).^2+(DX(:,:,2)-lambdawx(:,:,2)/beta2).^2+(DY(:,:,1)-lambdawy(:,:,1)/beta2).^2+(DY(:,:,2)-lambdawy(:,:,2)/beta2).^2);
s1(s1 == 0) = 1e-3;

s2 = sqrt((DX(:,:,3)-lambdawx(:,:,3)/beta2).^2+(DY(:,:,3)-lambdawy(:,:,3)/beta2).^2);
s2(s2 == 0) = 1e-3;

wx(:,:,3) = max(0, s2-eta*mu/beta2).*(DX(:,:,3)-lambdawx(:,:,3)/beta2)./s2;
wx(:,:,1) = max(0, s1-eta/beta2).*(DX(:,:,1)-lambdawx(:,:,1)/beta2)./s1;
wx(:,:,2) = max(0, s1-eta/beta2).*(DX(:,:,2)-lambdawx(:,:,2)/beta2)./s1;

wy(:,:,3) = max(0, s2-eta*mu/beta2).*(DY(:,:,3)-lambdawy(:,:,3)/beta2)./s2;
wy(:,:,1) = max(0, s1-eta/beta2).*(DY(:,:,1)-lambdawy(:,:,1)/beta2)./s1;
wy(:,:,2) = max(0, s1-eta/beta2).*(DY(:,:,2)-lambdawy(:,:,2)/beta2)./s1;
end

