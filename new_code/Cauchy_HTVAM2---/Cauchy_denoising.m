function [u,g1,g2] = Cauchy_denoising(x_n,alpha,M,r,beta,tol,maxiter,gamma,mu)
%% Description
%Solve the model£º
%argmin(u,g1,g2)=1/2£¨log(gamma^2+(u-x_n)^2)+mu*||u-u_0||.^2£©+alpha*(||g1-M/alpha||.^2+||g2-M/alpha||.^2)+g1.^2.*K(uTV)+g2.^2.*K(uLLT);
%x_n is the observation£¬x_t is the original image
%% 
g1=0;g2=0;u=x_n;psf=ones(r)/(r^2);
i=0;crit=1;
lamta.x=0;lamta.y=0;lamta.xx=0;lamta.xy=0;lamta.yx=0;lamta.yy=0;
p.x=0;p.y=0;q.xx=0;q.xy=0;q.yx=0;q.yy=0;lamta3=0;g=0;
u0=medfilt2(x_n);
%% 
while i < maxiter.out && crit > tol
    up=u;g1p=g1;g2p=g2;
    %% 
    [u,p,q,lamta]=ALM_nested(x_n,g,u,g1,g2,psf,beta,p,q,lamta,lamta3,tol,maxiter.in,gamma,mu,u0);
    %% 
    uTV=TotalVariation(u,'circular');
    uTV=imfilter(uTV,psf,'conv','circular');
    uLLT=LLT(u,'circular');
    uLLT=imfilter(uLLT,psf,'conv','circular');
    g1=M./(uTV+alpha);
    g2=M./(uLLT+alpha);
    %% 
    i=i+1;
%     u=min(max(0,u),1);
end

end

%% 
function [u,p,q,lamta]=ALM_nested(x_n,g,u,g1,g2,psf,beta,p,q,lamta,lamta3,tol,maxiter,gamma,mu,u0)
%% Description
%Solve the model£ºargmin(u)=1/2£¨log(gamma^2+(u-x_n)^2)+mu*||u-u_0||.^2£©+g1.^2.*K(uTV)+g2.^2.*K(uLLT)
%% 
kkmax =10;
kk=0;
i=0;crit=1;
sizeB=size(x_n);
otflaplace=psf2otf([0,1,0;1,-4,1;0,1,0],sizeB);
otfHstarH=psf2otf([0,0,1,0,0;0,2,-8,2,0;1,-8,20,-8,1;0,2,-8,2,0;0,0,1,0,0;],sizeB);
denom=1-otflaplace+otfHstarH;
g1temp=imfilter(g1.^2,psf,'conv','circular');
g2temp=imfilter(g2.^2,psf,'conv','circular');
%% 
while i < maxiter && crit > tol
    up=u;
    %% 
    nomin=-div(lamta.x+beta*p.x,lamta.y+beta*p.y,'circular')+...
        div2(lamta.xx+beta*q.xx,lamta.xy+beta*q.xy,lamta.yx+...
        beta*q.yx,lamta.yy+beta*q.yy,'circular')+beta*g+lamta3;
    nomin=fft2(nomin/beta);
    utemp=ifft2(nomin./denom);
    u=real(utemp);
    %% 
    [ux,uy]=grad(u,'circular');
    wx=ux-lamta.x/beta;
    wy=uy-lamta.y/beta;
    absw=sqrt(wx.^2+wy.^2);
    wtemp1=max(0,absw-g1temp/beta);
    absw(absw==0)=1;
    p.x=wtemp1.*wx./absw;
    p.y=wtemp1.*wy./absw;
    %%
    [uxx,uxy,uyx,uyy]=hessian(u,'circular');
    wxx=uxx-lamta.xx/beta;
    wxy=uxy-lamta.xy/beta;
    wyx=uyx-lamta.yx/beta;
    wyy=uyy-lamta.yy/beta;
    absw=sqrt(wxx.^2+wxy.^2+wyx.^2+wyy.^2);
    wtemp1=max(0,absw-g2temp/beta);
    absw(absw==0)=1;
    q.xx=wtemp1.*wxx./absw;
    q.xy=wtemp1.*wxy./absw;
    q.yx=wtemp1.*wyx./absw;
    q.yy=wtemp1.*wyy./absw;
    
    temp=gamma^2+(g-x_n).^2;
    fd =(g-x_n)./temp+mu*(g-u0)+beta*(g-u)+lamta3;
    while kk<kkmax && (norm(fd,'fro')>1e-5)
        fdd = (gamma^2-(g-x_n).^2)./(temp.^2)+mu+beta;
        g = g-fd./fdd;
        kk=kk+1;
        
        temp=gamma^2+(g-x_n).^2;
        fd =(g-x_n)./temp+mu*(g-u0)+beta*(g-u)+lamta3;       
    end
    kk=0;
    g(g<0)=0;
    %% 
    lamta.x=lamta.x+beta*(p.x-ux);
    lamta.y=lamta.y+beta*(p.y-uy);
    lamta.xx=lamta.xx+beta*(q.xx-uxx);
    lamta.xy=lamta.xy+beta*(q.xy-uxy);
    lamta.yx=lamta.yx+beta*(q.yx-uyx);
    lamta.yy=lamta.yy+beta*(q.yy-uyy);
    lamta3=lamta3+beta*(g-u);
    %% 
    i=i+1;
    crit=norm(u-up,'fro')/norm(up,'fro');
end
end

