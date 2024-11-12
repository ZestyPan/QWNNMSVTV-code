function [D]=TotalVariation(u,BoundaryCondition)
%������ȫ���, Ĭ��ѭ���߽�����
%ʾ����D=TotalVariation(u,'circular');
if nargin<2
    BoundaryCondition='circular';
end
[ux,uy]=grad(u,BoundaryCondition);
D=sqrt(ux.^2+uy.^2);
end