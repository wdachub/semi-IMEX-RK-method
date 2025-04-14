function L=  BC_CH_noflux_IMEX_L(L,D,D3,epsilon)

L([1,end],:)=D([1,end],:);%no flux for the scalar
L([2,end-1],:)=-epsilon*D3([1,end],:);%no flux for the chemical potential

end