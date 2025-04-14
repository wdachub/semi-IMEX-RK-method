function pt = SemiIMEXRKungeKutta2(odeEX, G, tspan, p0,M)
% SemiIMEXRKungeKutta2 performs second-order integration of a velocity field
% using semi implicit-explicit (IMEX) Runge-Kutta method.
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
if exist('M','var')
    existM=1;
else
    M=speye(length(p0));
    existM=0;

end

pt = p0;
for indt=1:length(tspan)-1
    pt=RKOneStep(odeEX,tspan(indt),dt, pt,G,M,existM);
end

end

function pt = RKOneStep(odeEX,t0,dt, p0,G,M,existM)
%
% GI=I-dt*G(p0);
% temp=odeEX(t0, p0);
% pt=GI\(p0+dt*temp);



GI=M-0.5*dt*G(p0);
f=odeEX(t0, p0);
if existM
    K1=GI\(M*p0+0.5*dt*f);
else
    K1=GI\(p0+0.5*dt*f);
end


GI=M-0.5*dt*G(K1);
f=odeEX(t0+0.5*dt, p0);
if existM
    K2=GI\(M*p0+0.5*dt*f);
else
    K2=GI\(p0+0.5*dt*f);
end
pt=2*K2-p0;


% L=I-0.5*dt*makeLinearOP(D,u0,testcase);
% K1=L\(u0+0.5*dt*f);
% L=I-0.5*dt*makeLinearOP(D,K1,testcase);
% K2=L\(u0+0.5*dt*f);
% u=2*K2-u0;


end

