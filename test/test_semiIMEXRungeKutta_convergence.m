clear
close all
addpath('..\tools')
addpath('..\RKmethod')

profile on
syms y(t)
eqn = diff(y,t) ==(-y+cos(t))*y+cos(t)*y;
% eqn = diff(y,t) ==(-y)*y;
cond = y(0) == 1;
f = dsolve(eqn,cond);
exact=matlabFunction(f, 'vars', t);
t0=0;
te=1;
odeEX=@(t,y)cos(t)*y;
odeIM=@(t,y)-y+cos(t);
% odeRK4=@(t,y)(-y+cos(t))*y+cos(t)*y;

% odeEX=@(t,y)myodeEX(t,y);
% odeIM=@(t,y)myodeIM(t,y);


% method=[2,2,1];
method=[2,3,4];
% method=[2,3,3];
order=method(1);


steplist=2.^(1:6);
dtlist=zeros(size(steplist));
ptlist=dtlist;
error=dtlist;
tic
for inds=1:length(steplist)

p0=ones(1);
nstep=steplist(inds);
tt=linspace(t0,te,nstep+1);
dtlist(inds)=tt(2)-tt(1);
pt=semiIMEXRungeKutta(odeEX,odeIM, tt, p0, method,0);

ptlist(inds)=pt; 
error(inds)=max(abs(pt-exact(te)),[],'all')/max(abs(exact(te)));
% pt
end
toc
error(end)

 profile viewer

figure(1)
movegui('northwest')
loglog(dtlist,error,'r')
hold on
loglog(dtlist,error(end)/dtlist(end)^order*dtlist.^order,'--b')
hold off
title('convergence rate')
legend('numerical convergence rate',['order=',num2str(order)],'Location','southeast')

