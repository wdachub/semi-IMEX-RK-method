function pt = semiIMEXRungeKuttaTrapezoid(varargin)
% semiIMEXRungeKuttaTrapezoid Performs time integration using the trapezoid method, a Semi-Implicit-Explicit (IMEX) Runge-Kutta method.
%
% This function integrates a velocity field governed by equations of the form:
%   y' = f(t, y) + G(t, y)y,
% where y is a vector, f(t, y) is a vector-valued function (explicit part), and G(t, y) is a matrix (implicit part).
% The explicit term f(t, y) is treated with an explicit method, while the implicit term G(t, y)y is handled implicitly.
% Trapezoid method is used.
%
% Inputs:
%   odeEX  : Function handle for the explicit ODE part, returning f(t, y).
%   odeIM  : Function handle for the implicit ODE part, returning G(t, y).
%   tspan  : Vector specifying the time interval for integration [t0, te].
%   pt     : Initial position of the point, which can be in any orientation.
%   method : Vector of three components specifying the numerical method:
%            - method(1): Order of the method.
%            - method(2): Number of stages.
%            - method(3): Index (used when multiple methods with the same order and stage count exist).
%   Gtype  : Scalar indicating the dependence of G(t, y):
%            - 0: G(t, y) is independent of time (i.e., G(y)).
%            - 1: G(t, y) depends on time.
%
% Outputs:
%   pt     : Final position of the point after integration.
%
%
% Butcher Tableau Structure:
% - BT.s: Number of stages.
% - BT.aE, BT.bE, BT.cE: Coefficients for the explicit method.
% - BT.aI, BT.bI, BT.cI: Coefficients for the implicit method.
% - BT.alpha: Special parameter for handling solution reconstruction:
%   - If BT.alpha = NaN, the solution is computed as a linear combination of intermediate steps.
%   - If BT.alpha = 1, the solution corresponds to the last stage.
%   - If BT.alpha is finite, the solution is computed using the last stage value.

%% Zero time interval check
if abs(varargin{3}(end)-varargin{3}(1))<eps
    return
end


%% Input validation
p = inputParser;
addRequired(p, 'odeEX', @(x) isa(x, 'function_handle'));
addRequired(p, 'odeIM', @(x) isa(x, 'function_handle'));
addRequired(p, 'tspan', @(x) isvector(x));
addRequired(p, 'pt', @isnumeric);


% Parse inputs
parse(p, varargin{:});

% Extract values from parsed input
odeEX = p.Results.odeEX;
odeIM = p.Results.odeIM;
tspan = p.Results.tspan;
pt = p.Results.pt;

%% Butcher tableau for Trapezoid 
% 
% %A stable
% BT.alpha=0.5;
% BT.s=3;%stage
% BT.aE = zeros(BT.s);
% BT.aE(2,1) = 0.5;
% BT.aE(3,2) = 0.5;
% BT.bE =[0,1,0];
% BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
% %implicit
% BT.aI = zeros(BT.s);
% BT.aI(2,2) = 0.5;
% BT.aI(3,3) = 0.5;
% BT.bI =[0,0,0,1];
% BT.cI= sum(BT.aI,2);%the method satisfy row sum condition


%%
if any(size(pt)>10)
    %check if the input length is large, if so, use sparse matrix
    %I am not sure if the threshold 10 is optimal.
    I=speye(length(pt));
else
    I=eye(length(pt));
end



%% Time integration loop
for indt=1:length(tspan)-1
    dt=tspan(indt+1)-tspan(indt);
    if abs(dt)>eps
        pt=RKOneStep(odeEX,odeIM, tspan(indt),dt, pt,I);
    end
end

end
function pt = RKOneStep(odeEX,odeIM, t0,dt, pt,I)
% RKOneStep performs one step of the Trapezoid method, a
% Semi-Implicit-Explicit (IMEX) Runge-Kutta method.
% Inputs:
%   - odeEX: External function for the explicit part of the ODE (f(y))
%   - odeIM: Implicit function for the implicit part of the ODE (G)
%   - t0: Current time
%   - dt: Time step
%   - p0: Initial condition at time t0
%   - BT: A structure containing coefficients for the Runge-Kutta method
%
% Output:
%   - pt: Solution at the next time step (t0 + dt)



L=I-dt*0.5*odeIM(t0+0.5*dt,pt);
K2=L\(pt+0.5*dt*odeEX(t0,pt));%% Solve the equation L * K{indi} = K{indi}


L=I-dt*0.5*odeIM(t0+0.5*dt,K2);
K3=L\(pt+0.5*dt*odeEX(t0+0.5*dt,K2));%% Solve the equation L * K{indi} = K{indi}
pt=2*K3-pt;


end