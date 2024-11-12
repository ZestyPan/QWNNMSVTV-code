function [u]=div2(p11,p12,p21,p22,BoundaryCondition)
%用于求Hessian算子的共轭算子, 默认循环边界
%使用循环边界时函数使用diff故运算速度快, 其他边界条件使用imfilter
%u=D_x^+(D_x^-(p11))+D_y^-(D_x^-(p12))+D_x^-(D_y^-(p21))+D_y^+(D_y^-(p22));
%示例: u=div2(p11,p12,p21,p22);
if nargin<5
    BoundaryCondition='circular';
end
if strcmp(BoundaryCondition,'circular')
    u=Dxforward(Dxbackward(p11))+Dybackward(Dxbackward(p12))+...
        Dxbackward(Dybackward(p21))+Dyforward(Dybackward(p22));
else
    u11=imfilter(p11,[1,-2,1],BoundaryCondition);
    u12=imfilter(p12,[1 -1 0;-1 1 0;0 0 0],BoundaryCondition);
    u21=imfilter(p21,[1 -1 0;-1 1 0;0 0 0],BoundaryCondition);
    u22=imfilter(p22,[1;-2;1],BoundaryCondition);
    u=u11+u12+u21+u22;
end
end