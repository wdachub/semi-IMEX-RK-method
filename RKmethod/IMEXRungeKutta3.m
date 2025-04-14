function pt = IMEXRungeKutta3(odeEX, odeIM, tspan, p0,BC)
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
%explicit
ExRK.a = zeros(4,4);
ExRK.a(3,2)=1;
ExRK.a(4,2:3)=[1/4,1/4];
ExRK.b=[0,1/6,1/6,2/3];
ExRK.c=[0,0,1,1/2];
ExRK.s=3;%stage
%implicit
%     alpha=0.24169426078821;%from paper, only have 14 digits.seems not satisfy L-stable condition
%     beta=0.06042356519705;
%     eta=0.12915286960590;
%     beta=(27-sqrt(409))/320;%derived by myself, wrong,I solve wrong equation for L-stable.
%     alpha=4*beta;
%     eta=(13+sqrt(409))/160;

%     beta=(27+sqrt(409))/320;
%     alpha=4*beta;
%     eta=(13-sqrt(409))/160;%eta is nagetive


%     beta=(9-sqrt(57))/24;
%     alpha=4*beta;
%     eta=(-6+sqrt(57))/12;

beta=1/8;
alpha=4*beta;
eta=0;

ImRK.a= zeros(4,4);
ImRK.a(1,1)=alpha;
ImRK.a(2,1:2)=[-alpha,alpha];
ImRK.a(3,2:3)=[1-alpha,alpha];
ImRK.a(4,1:4)=[beta,eta,1/2-beta-eta-alpha,alpha];
ImRK.b=[0,1/6,1/6,2/3];
ImRK.c=[alpha,0,1,1/2];
ImRK.s=4;

%%


pt = p0;

G=odeIM;
L=1-dt*ImRK.a(1,1)*G;


for indt=1:length(tspan)-1
    pt=RKOneStep(odeEX, tspan(indt),dt, pt,ExRK,ImRK,G,L,BC);
end

end

function p0 = RKOneStep(odeEX,t0,dt, p0,ExRK,ImRK,G,L,BC)
maxs=max(ExRK.s,ImRK.s);

[Us{1:maxs}] = deal(p0);
for indj=1:maxs


    if exist('BC','var')
        [L,Us{indj}] = BC(L,Us{indj});

    end
    Utemp=L\Us{indj};%solve equation

    if indj>1
        Ftemp=odeEX(t0+ExRK.c(indj)*dt, Utemp);%EX is 3 stage, while IM is 4 stage.
    end
    Gtemp=G*Utemp;%replace G  with g(y)

    for indi=indj+1:maxs
        if abs(ExRK.a(indi,indj))>eps
            Us{indi} = Us{indi} + dt*ExRK.a(indi,indj)*Ftemp;
        end
        if abs(ImRK.a(indi,indj))>eps
            Us{indi} =Us{indi} + dt*ImRK.a(indi,indj)*Gtemp;
        end
    end


    %     if abs(ExRK.b(indj))>eps%ExRK.b==ImRK.b
    %         p0 = p0 + dt*ExRK.b(indj)*(Ftemp+Gtemp);
    %     end

    if abs(ExRK.b(indj))>eps
        p0 = p0 + dt*ExRK.b(indj)*Ftemp;
    end
    if abs(ImRK.b(indj))>eps
        p0 = p0 + dt*ImRK.b(indj)*Gtemp;
    end
end

end

