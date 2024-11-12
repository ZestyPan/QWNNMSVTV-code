function [D]=LLT(u,BoundaryCondition)
%������߽�ȫ��ʹ��LLTģ��, Ĭ��ѭ���߽�����
%���ж���ƫ�����õ����ĸ�ʽ�����ƫ�����õ���ǰ��ǰ��ʽ
%ʾ����D=LLT(u,'circular');
if nargin<2
    BoundaryCondition='circular';
end
[pxx,pxy,pyx,pyy]=hessian(u,BoundaryCondition);
D=sqrt(pxx.^2+pxy.^2+pyx.^2+pxx.^2+pyy.^2);
end