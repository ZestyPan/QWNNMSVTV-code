function [D]=LLT(u,BoundaryCondition)
%用于求高阶全变差，使用LLT模型, 默认循环边界条件
%其中二阶偏导数用的中心格式，混合偏导数用的向前向前格式
%示例：D=LLT(u,'circular');
if nargin<2
    BoundaryCondition='circular';
end
[pxx,pxy,pyx,pyy]=hessian(u,BoundaryCondition);
D=sqrt(pxx.^2+pxy.^2+pyx.^2+pxx.^2+pyy.^2);
end