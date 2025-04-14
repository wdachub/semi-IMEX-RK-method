%If the exact solution is not available, we compute a reference solution using a high-order method with a fine time step size and treat it as the exact solution for error estimation.

clear
close all

addpath('..\tools')
addpath('..\RKmethod')

M=2^7;
Lx=2*pi;
x=linspace(-Lx/2,Lx/2,M+1)';
ns=5;
D=diffmat(x,ns,1);
D=D{1};
I=speye(length(x));
t0=0;
te=1;


odeEX=@(t,y)cos(x)*sin(t);
odeIM=@(t,y)D*spdiags(1+y.^2,0,length(y),length(y))*D;
cIC=@(x,t0)zeros(size(x));


%%
%
methodlist={[1,2,1],[2,3,1],[2,3,4],[3,5,1],[3,5,2]};
steplist=2.^(4:7);


ct=cIC(x,t0);
method=[3,5,2];
tt=linspace(t0,te,steplist(end)*4+1);
ce=semiIMEXRungeKutta(odeEX,odeIM, tt, ct, method,0,@(L,rhs,pt)BC_nonlinear_diffusion(L,rhs,pt,D));
cemax=max(abs(ce));
plot(x,ce,'-r')
% return
errorlist=zeros(length(methodlist),length(steplist)-1);


for indm=1:length(methodlist)
    method=methodlist{indm};
    order=method(1);

    dtlist=zeros(size(steplist));
    tic
    for inds=1:length(steplist)
        tic
        ct=cIC(x,t0);

        nstep=steplist(inds);


        tt=linspace(t0,te,nstep+1);
        dtlist(inds)=tt(2)-tt(1);
        ct=semiIMEXRungeKutta(odeEX,odeIM, tt, ct, method,0,@(L,rhs,pt)BC_nonlinear_diffusion(L,rhs,pt,D));
        errorlist(indm,inds)=max(abs(ct-ce))/cemax;

        toc
    end

end
to_latex_convergence_table(steplist, errorlist, 'test.txt')

