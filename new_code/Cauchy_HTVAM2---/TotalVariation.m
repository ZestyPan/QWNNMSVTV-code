function [D]=TotalVariation(u,BoundaryCondition)
%用于求全变差, 默认循环边界条件
%示例：D=TotalVariation(u,'circular');
if nargin<2
    BoundaryCondition='circular';
end
[ux,uy]=grad(u,BoundaryCondition);
D=sqrt(ux.^2+uy.^2);
end