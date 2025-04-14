function pt= RungeKutta4(odefun, tspan, p0)
% Runge-Kutta integrating velocity field.
% using initial guess
% Input:
% p0     : the initial position of a point.  vertical or horizontal
% velocity field is important
% [t0,te]: the time interval, assuming uniform time step size
% vel    : the velocity history,
% order  : order of the integrator
%
% Output:
% pt     : the ending position of the point

dt = tspan(2)-tspan(1);
if abs(dt)<eps
    pt = p0;
    return
end
pt = p0;

for indt=1:length(tspan)-1
    pt=RKOneStep(odefun, tspan(indt),dt, pt);

end

end

function pt = RKOneStep(odefun, t0,dt, p0)
k1=odefun(t0, p0);
k2=odefun(t0+0.5*dt, p0+dt*0.5*k1);
k3=odefun(t0+0.5*dt, p0+dt*0.5*k2);
k4=odefun(t0+dt, p0+dt*k3);
pt=p0+dt*(k1+2*k2+2*k3+k4)/6;
end