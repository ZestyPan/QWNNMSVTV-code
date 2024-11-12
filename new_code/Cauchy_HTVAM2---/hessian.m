function [uxx,uxy,uyx,uyy]=hessian(u,BoundaryCondition)
%用于求Hessian矩阵, 默认循环边界
%使用循环边界时函数使用diff故运算速度快, 其他边界条件使用imfilter
%uxx=D_x^-(D_x^+(u));
%uxy=D_x^+(D_y^+(u));
%uyx=D_y^+(D_x^+(u));
%uyy=D_y^-(D_y^+(u));
%示例：[uxx,uxy,uyx,uyy]=hessian(u);
if nargin<2
    BoundaryCondition='circular';
end
if strcmp(BoundaryCondition,'circular')
    uxx=Dxbackward(Dxforward(u));
    uxy=Dxforward(Dyforward(u));
    uyx=uxy;
    uyy=Dybackward(Dyforward(u));
else
    uxx=imfilter(u,[1,-2,1],BoundaryCondition);
    uxy=imfilter(u,[0 0 0;0 1 -1;0 -1 1],BoundaryCondition);
    uyx=uxy;
    uyy=imfilter(u,[1;-2;1],BoundaryCondition);
end
end