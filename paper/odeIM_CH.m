function dy = odeIM_CH(t,y,D,D4)
dy=-epsilon^2*D4+D*spdiags(3*y.^2-1,0,length(y),length(y))*D;

dy

end