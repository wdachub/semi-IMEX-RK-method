function [L,rhs]=  BC_nonlinear_diffusion(L,rhs,pt,D)

% rhs([1,2,end-1,end])=[0,0,0,0];
rhs([1,end])=0;
%periodic boundary
L(1,:)=0;
L(1,1)=1;
L(1,end)=-1;% y(1)=y(end);

%dy=D*y, dy(1)=dy(end)
L(end,:)=D(1,:)-D(end,:);

%mean zero
% rhs(2)=0;
% L(2,:)=1;
end