function pt = IMEXRungeKutta2(odeEX, odeIM, tspan, pt,BCL,BCRhs)
% IMEXRungeKutta3 performs third-order integration of a velocity field
% using implicit-explicit (IMEX) Runge-Kutta method.
%
% Inputs:
%   odeEX  : function handle for the explicit ODE part.
%   G      : function handle for the implicit ODE part.
%   tspan  : vector specifying the time interval for integration [t0, te].
%   p0     : initial position of a point, can be in any orientation.
%
% Outputs:
%   pt     : final position of the point after integration.
%
% Example of use:
%   pt = IMEXRungeKutta3(@odeFuncEX, @odeFuncIM, [0, 10], [1; 0]);
%
% This function integrates a velocity field that is specified by the user.
% The integration uses a third-order IMEX Runge-Kutta method, which handles
% the dynamics of the system through both implicit and explicit stages.


dt = tspan(2)-tspan(1);
if abs(dt)<eps
    pt = p0;
    return
end

Nt=ceil((tspan(end)-tspan(1))/dt);
tspan=linspace(tspan(1),tspan(end),Nt+1);


%%

%%

gamma=(2-sqrt(2))/2;
delta=-1/sqrt(2);

BT.alpha=1;
BT.s=3;%stage
BT.aE = zeros(BT.s);
BT.aE(2,1) = gamma;
BT.aE(3,[1,2]) = [delta, 1-delta];
BT.bE =[delta, 1-delta,0];
BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
%implicit
BT.aI = zeros(BT.s);
BT.aI(2,2) = gamma;
BT.aI(3,[2,3]) = [1-gamma,gamma];
BT.bI =[0,1-gamma,gamma];
BT.cI= sum(BT.aI,2);%the method satisfy row sum condition



G=odeIM;%assume odeIM is time-independent.
L= speye(length(pt))-dt*gamma*G;
if exist('BCL','var')
    L= BCL(L);
end

for indt=1:length(tspan)-1
    pt=RKOneStep(odeEX, tspan(indt),dt, pt,G,L,BT,BCRhs);
    if any(isnan(pt))
        return
    end

end

end

function pt = RKOneStep(odeEX,t0,dt, p0,G,L,BT,BCRhs)

Fs{1}=dt*odeEX(t0,p0);

K{1}=p0+BT.aE(2,1)*Fs{1};

K{1}=L\BCRhs(K{1});
Fs{2}=dt*odeEX(t0+BT.cE(2)*dt,K{1});
Gs{1}=dt*G*K{1};
K{2}=p0+BT.aE(3,1)*Fs{1}+BT.aE(3,2)*Fs{2}+ BT.aI(3,2)*Gs{1};
pt=L\BCRhs(K{2});%pt=K{2};



end

