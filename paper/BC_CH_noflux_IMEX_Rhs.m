function rhs=  BC_CH_noflux_IMEX_Rhs(rhs,ct,D,fx)


rhs([1,end])=0;

temp=-D*fx(ct);

rhs([2,end-1],:)=temp([1,end]);%no flux for the chemical potential

end