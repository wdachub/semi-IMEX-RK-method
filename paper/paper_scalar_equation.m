clear
close all
addpath('..\tools')
addpath('..\RKmethod')

syms y(t)
eqn = diff(y,t) ==(-y+cos(t))*y+cos(t)*y;
cond = y(0) == 1;
f = dsolve(eqn,cond);
exact=matlabFunction(f, 'vars', t);
t0=0;
te=0.5;
odeEX=@(t,y)cos(t)*y;
odeIM=@(t,y)-y+cos(t);
% odeEX=@(t,y)myodeEX(t,y);
% odeIM=@(t,y)myodeIM(t,y);


methodlist={[3,5,1]};
steplist=2.^(2:6);


errorlist=zeros(length(methodlist),length(steplist));
for indm=1:length(methodlist)
    method=methodlist{indm};
    order=method(1);

    dtlist=zeros(size(steplist));
    ptlist=dtlist;
    tic
    for inds=1:length(steplist)
        p0=ones(1);
        nstep=steplist(inds);
        tt=linspace(t0,te,nstep+1);
        dtlist(inds)=tt(2)-tt(1);
        pt=semiIMEXRungeKutta(odeEX,odeIM, tt, p0, method);
        ptlist(inds)=pt;
        errorlist(indm,inds)=max(abs(pt-exact(te)),[],'all')/max(abs(exact(te)));
        
    end
    toc

end
%save result to txt file. 
to_latex_convergence_table(steplist, errorlist, 'test.txt')


