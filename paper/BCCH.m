function [L,rhs]=  BCCH(L,rhs,pt,D,D3,epsilon)
rhs([1,2,end-1,end])=[0,0,0,0];
L([1,end],:)=D([1,end],:);%no flux for the scalar
temp=-epsilon^2*D3+spdiags(3*pt.^2-1,0,length(pt),length(pt))*D;
L([2,end-1],:)=temp([1,end],:);%no flux for the chemical potential
end