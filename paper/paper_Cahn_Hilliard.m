%If the exact solution is not available, we compute a reference solution using a high-order method with a fine time step size and treat it as the exact solution for error estimation.

clear
close all

addpath('..\tools')
addpath('..\RKmethod')

M=2^7;
Lx=40;
x=linspace(-1,1,M)';
x=sign(x).*abs(x).^1.5;
x=Lx/2*x;
ns=5;

D=diffmat(x,ns,4);
cIC=@(x,t0)tanh(x);

D4=D{4};
D3=D{3};
D2=D{2};
D=D{1};
I=speye(length(x));
t0=0;
te=1;

epsilon=1;
odeEX=@(t,y)zeros(size(y));
odeIM=@(t,y)-epsilon^2*D4+D*spdiags(3*y.^2-1,0,length(y),length(y))*D;



%%
% methodlist={[2,2,1],[2,3,1],[2,3,2],[2,3,3],[3,4,1],[3,5,1],[2,2,1]};
methodlist={[1,2,1],[2,3,1],[2,3,4],[3,5,1],[3,5,2]};
steplist=2.^(8:11);

ct=cIC(x,t0);
method=[3,5,2];
tt=linspace(t0,te,steplist(end)*4+1);
ce=semiIMEXRungeKutta(odeEX,odeIM, tt, ct, method,0,@(L,rhs,pt)BCCH(L,rhs,pt,D,D3,epsilon));
cemax=max(abs(ce));



errorlist=zeros(length(methodlist),length(steplist));


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
        ct=semiIMEXRungeKutta(odeEX,odeIM, tt, ct, method,0,@(L,rhs,pt)BCCH(L,rhs,pt,D,D3,epsilon));
        errorlist(indm,inds)=max(abs(ct-ce))/cemax;

        toc
    end

end
to_latex_convergence_table(steplist, errorlist, 'test.txt')

