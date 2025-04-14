function pt = semiIMEXRungeKutta(varargin)
% SEMIIMEXRUNGEKUTTA Performs time integration using a Semi-Implicit-Explicit (IMEX) Runge-Kutta method.
%
% This function integrates a velocity field governed by equations of the form:
%   y' = f(t, y) + G(t, y)y,
% where y is a vector, f(t, y) is a vector-valued function (explicit part), and G(t, y) is a matrix (implicit part).
% The explicit term f(t, y) is treated with an explicit method, while the implicit term G(t, y)y is handled implicitly.
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
% Notes:
% - Designed to handle different methods with given Butcher tableaux.
% - Uses multiple conditional statements to skip zero elements in the Butcher tableau, making it more general but slower for scalar equations compared to methods tailored for specific method.
%
% Butcher Tableau Structure:
% - BT.s: Number of stages.
% - BT.aE, BT.bE, BT.cE: Coefficients for the explicit method.
% - BT.aI, BT.bI, BT.cI: Coefficients for the implicit method.
% - BT.alpha: Special parameter for handling solution reconstruction:
%   - If BT.alpha = NaN, the solution is computed as a linear combination of intermediate steps.
%   - If BT.alpha = 1, the solution corresponds to the last stage.
%   - If BT.alpha is finite, the solution is computed using the last stage value.
% 
% I assume for each row, there is at least one zero element. in
%  other word, K{indi} must be different from K{indi-1}. A table voilate
%  this property will lead to additional calculation
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

% Default values for optional parameters
addOptional(p, 'method', [1, 2, 1], @isvector); % Default method: Forward-backward Euler
addOptional(p, 'Gtype', 1, @isnumeric); % Default: time-dependent G
addOptional(p, 'BC', [], @(x) isa(x, 'function_handle'));

% Parse inputs
parse(p, varargin{:});



% Extract values from parsed input
odeEX = p.Results.odeEX;
odeIM = p.Results.odeIM;
tspan = p.Results.tspan;
pt = p.Results.pt;
method = p.Results.method;
Gtype = p.Results.Gtype;

% Check if 'BC' was provided
if isempty(p.Results.BC)
    BT.BCExists = false;
     BC = [];
else
    BT.BCExists = true;
    BC = p.Results.BC;
end



if length(method) < 3
    method(3) = 1;
end

order = method(1);
stage = method(2);
index = method(3);

BT.Gtype = Gtype; % Time-dependent flag for G

%% Butcher tableau initialization
if order==1
    % First-order method, Forward-backward Euler, L stable
    %explicit
    BT.alpha=1;
    BT.s=2;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 1;
    BT.bE =[1,0];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,2) = 1;
    BT.bI =[0,0,1];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition

elseif order==2&&stage==2
    %A stable
    BT.s=2;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 0.5;
    BT.bE =[0,1];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,2) = 0.5;
    BT.bI =[0,1,0];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
elseif order==2&&stage==3&&index==1
    %A stable
    BT.alpha=0.5;
    BT.s=3;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 0.5;
    BT.aE(3,2) = 0.5;
    BT.bE =[0,1,0];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,2) = 0.5;
    BT.aI(3,3) = 0.5;
    BT.bI =[0,0,0,1];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
elseif order==2&&stage==3&&index==2
    %A stable
    BT.alpha=1;
    BT.s=3;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 1;
    BT.aE(3,[1,2]) = [0.5,0.5];
    BT.bE =[0.5,0.5,0];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,2) = 1;
    BT.aI(3,[1,3]) = [0.5,0.5];
    BT.bI =[0.5,0,0,0.5];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
elseif order==2&&stage==3&&index==3
    %L stable
    gamma=1-1/sqrt(2);
    BT.s=3;%stage
    BT.aE = zeros(BT.s);
    BT.aE(3,1) = 1;
    BT.bE =[0.5,0,0.5];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(1,1) = gamma;
    BT.aI(2,1) = 1-gamma;
    BT.aI(3,[1,3]) = [1-2*gamma,gamma];
    BT.bI =[0.5,0,0.5,0];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
elseif order==2&&stage==3&&index==4
    %L stable
     BT.alpha=1;
    gamma=1-1/sqrt(2);
    BT.s=3;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 1;
    BT.aE(3,[1,2]) = [0.5,0.5];
    BT.bE =[0.5,0.5,0];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    
    BT.aI(2,[1,2]) = [1-gamma,gamma];
    BT.aI(3,[1,2,3]) = [0.5,0.5-gamma,gamma];
    BT.bI =[0.5,0.5-gamma,0,gamma];
    % BT.bI =[0.5,0,0,0.5];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
elseif order==3&&stage==4&&index==1
    %L stable  since cE(3)=cE(4), it reduceds one evaluation of G
    BT.s=4;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 0.7775079538595848;
    BT.aE(3,[1,2]) = [0.3850382624054263,0.2733484980719337];
    BT.aE(4,[1,2,3]) = [0.2905474198112961,0.1784065415104640,0.1894327991556034];
    BT.bE =[0.2486553715043413,0.04469938464765911, 0.3828282521031255,0.3238169917448679];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,1:2) = [0.5668275181562270,0.2106804357033578];
    BT.aI(3,[1,2,3]) = [0.3481097445529071,0.1497169356151823,0.1605600803092672];
    BT.aI(4,:) = [0.3299758037920577,0.1113697479208660,0.1255619659848192,0.09147924277961349];
    BT.bI =[0.2486553715043413,0.04469938464765911, 0.3828282521031255,0.3238169917448679,0];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
elseif order==3&&stage==4&&index==2
    %L stable
    BT.s=4;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 0.1674644883689973;
    BT.aE(3,[1,2]) = [0.4653041332000116,0.4677803889102060];
    BT.aE(4,[1,2,3]) = [0.2071239388016568, 0.2087440870871372, 0.2658255929491791];
    BT.bE =[0.08287386878034993, 0.2776394645104248,...
        0.06989620096971407, 0.5695904657395114];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,2) = 0.1674644883689973;
    BT.aI(3,[1,2,3]) = [0.3228096995961422,0.4274100950975858,0.1828647274164901];
    BT.aI(4,:) = [0.2560291354472103,0.1529665372918773,0.1499795496788363,0.1227183964200491];
    BT.bI =[0.08287386878034993, 0.2776394645104248,...
        0.06989620096971407, 0.5695904657395114,0];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition

elseif order==3&&stage==4&&index==3
    %not A stable, but has less function envaluation
    BT.s=4;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 0.1704760463630139;
    BT.aE(3,[1,2]) = [0.2970657136981410,0.5089374297855375];
    BT.aE(4,[1,2,3]) = [0.07234929368319974,0.08738793876247758,0.4197565976551324];
    BT.bE =[0.1410668830792905,0.1606899298515207,0.3001139168928477,0.3981292701763415];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,2) = 0.1704760463630139;
    BT.aI(3,[2,3]) = [0.7237525628681173,0.08225058061556137];
    BT.aI(4,2:4) = [0.2392125577156714,0.1145362811310622,0.2257449912540760];
    BT.bI =[0.1410668830792905,0.1606899298515207,0.3001139168928477,0.3981292701763415,0];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
elseif order==3&&stage==5&&index==1
    %L stable 3 inversion
    BT.alpha=1;
    BT.s=5;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 0.64116921315526898;
    BT.aE(3,[1,2]) = [0.39058950600403958,0.86314276923850819];
    BT.aE(4,[1,2,3]) = [0.42747115807408170,0.35555178088542744,0.21697706104049089];
    BT.aE(5,1:4) =[0.30991530721474964,0.32596239153256790,-0.28817520861282836,0.65229750986551083];
    BT.bE =[0.30991530721474964,0.32596239153256790,-0.28817520861282836,0.65229750986551083,0];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,1:2) = [0.30312000893712265,0.33804920421814655];
    BT.aI(3,1:3) = [0.39058950600403963, 0.46290999159550344, 0.40023277764300441];
    BT.aI(4,1:3) = [0.43415392037526129,0.34187417721762819,0.22397190240711046];
    BT.aI(5,1:5) = [0.30991530721474964,0.32596239153256790,-0.28817520861282836,0,0.65229750986551083];
    BT.bI =[0.30991530721474964,0.32596239153256790,-0.28817520861282836,0.65229750986551083,0,0.65229750986551083];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
elseif order==3&&stage==5&&index==2
     %L stable 4 inversion
    BT.alpha=1;
    BT.s=5;%stage
    BT.aE = zeros(BT.s);
    BT.aE(2,1) = 0.37729778462711194;
    BT.aE(3,[1,2]) = [0.32109244734547510,0.67890755265452751];
    BT.aE(4,[1,2,3]) = [0.29583591899535783,0.32786792139864995,0.37629615960599228];
    BT.aE(5,1:4) =[0.058262270658744675,0.70938840176878493,-0.20706199805500403,0.43941132562747443];
    BT.bE =[0.058262270658744675,0.70938840176878493,-0.20706199805500403,0.43941132562747443,0];
    BT.cE= sum(BT.aE,2);%the method satisfy row sum condition
    %implicit
    BT.aI = zeros(BT.s);
    BT.aI(2,1:2) = [0.27090231391056940,0.10639547071654235];
    BT.aI(3,1:3) = [0.32109244734547354, 0.45805080731378267, 0.22085674534074654];
    BT.aI(4,1:4) = [0.44587480986461181,0.086919861210029870,0.33728474074652454,0.12992058817883403];
    BT.aI(5,1:5) = [0.058262270658745036,0.70938840176878437,-0.20706199805500353,-0.21780858432897851,0.65721990995645263];
    BT.bI =[0.058262270658745036,0.70938840176878437,-0.20706199805500353,-0.21780858432897851,0,0.65721990995645263];
    BT.cI= sum(BT.aI,2);%the method satisfy row sum condition
else

    error('Unsupported method configuration');

end

%% Determine method type (alpha value)
if ~isfield(BT, 'alpha')
    if BT.bE(end)==0

        temp=[BT.aE(BT.s,1:end-1),BT.aI(BT.s,1:end)]./[BT.bE(1:end-1),BT.bI([1:end-2,end])];
        temp = temp(isfinite(temp)); % Keeps only finite values

        if  max(abs(temp-temp(1))) < 10*eps
            %check if all element are the same
            BT.alpha=temp(1);
        else
            BT.alpha=nan;
        end

    else

        BT.alpha=nan;
    end
    clear temp;

end
% S = rmfield(S, 'fieldname');
BT.is_need_update=false(1,BT.s);
BT.is_need_Gs=false(1,BT.s);
BT.is_need_Fs=false(1,BT.s);
BT.is_need_RGs=false(1,BT.s);
BT.is_need_solve=false(1,BT.s);
for indi=1:BT.s
    BT.is_need_solve(indi)=logical(BT.aI(indi,indi));
    % Check if there are any updates for the current stage (indi)
    BT.aEp{indi}= find(BT.aE(indi,1:indi-1));  % Stages with nonzero entries
    BT.aIp{indi}= find(BT.aI(indi,1:indi-1));  % Stages with nonzero entries
  
    BT.is_need_update(indi)=any(BT.aEp{indi})||any(BT.aIp{indi});
    % Check if we need to compute Gs for the current stage
    if indi==BT.s
        BT.is_need_Gs(indi)=any(BT.bI(end-1:end))&& isnan(BT.alpha);
    else
        %determine if Gs will be used in computing Ks by checking the element in aI.
        %determine if Gs will be used in computing pt when alpha is notfinite.
        BT.is_need_Gs(indi)=any(BT.aI(indi+1:end,indi))|| (isnan(BT.alpha)&&logical(BT.bI(indi)));
    end
    if indi>1
        BT.is_need_RGs(indi)=(~BT.is_need_Gs(indi-1)) || (BT.Gtype==1 && abs(BT.cI(indi)-BT.cI(indi-1))>20*eps);
    end
    BT.is_need_Fs(indi)=any(BT.aE(indi+1:end,indi))|| (isnan(BT.alpha)&&logical(BT.bE(indi)));% Check if we need to compute the explicit function for future stages
end

%%
if any(size(pt)>10)
    %check if the input length is large, if so, use sparse matrix
    %I am not sure if the threshold 10 is optimal.
    BT.I=speye(length(pt));
else
    BT.I=eye(length(pt));
end



%% Time integration loop
for indt=1:length(tspan)-1
    dt=tspan(indt+1)-tspan(indt);
    if abs(dt)>eps
        pt=RKOneStep(odeEX,odeIM, tspan(indt),dt, pt,BT,BC);
    end
end

end
function pt = RKOneStep(odeEX,odeIM, t0,dt, p0,BT,BC)
% RKOneStep performs one step of a Runge-Kutta method for an ODE system.
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

sz = size(p0); % Get the size of the initial condition vector
pt = zeros(sz); % Initialize the solution vector for the next time step
K = cell(1, BT.s); % Cell array to store intermediate results for K
Fs = cell(1, BT.s); % Cell array to store explicit function evaluations
Gs = cell(1, BT.s); % Cell array to store implicit function evaluations

for indi=1:BT.s

    if BT.is_need_update(indi)    % Check if there are any updates for the current stage (indi)
        Ktemp = zeros(sz);% Temporary variable to accumulate K values
        for indj=BT.aEp{indi}
            Ktemp = Ktemp + BT.aE(indi,indj)*Fs{indj};% Update with explicit term
        end
        for indj=BT.aIp{indi}
            Ktemp = Ktemp + BT.aI(indi,indj)*Gs{indj}; % Update with implicit term
        end
        K{indi}=p0+dt*Ktemp;% Update K for current stage
    else
        K{indi}=p0; % If no updates, set K to initial condition
    end

    % Determine whether we need to solve the implicit system for this stage
    if BT.is_need_solve(indi)
        % Compute matrix G. It is evaluate at K{indi-1}. K{0}=p0
        if indi==1
            GM=odeIM(t0+dt*BT.cI(indi),p0);
        else
            % Reuse the previous GM from the previous stage if the condition holds
            if  BT.is_need_RGs(indi)
                % compute GM if GM doesn't exist, OR,
                % if GM is time depedent and it is evaluated at different t
                % in previous stage

                GM=odeIM(t0+dt*BT.cI(indi),K{indi-1});%recompute G if the condition doesn't hold
            end
        end

        L=BT.I-dt*BT.aI(indi,indi)*GM;
        if  BT.BCExists
            if indi==1
                [L,K{indi}]=BC(L,K{indi},p0);%add boundary condition
            else
                [L,K{indi}]=BC(L,K{indi},K{indi-1});%add boundary condition

            end

        end
        K{indi}=L\K{indi};%% Solve the equation L * K{indi} = K{indi}
    end


    if BT.is_need_Gs(indi)
        %meaning GM has been computed

        if indi==BT.s
            if logical(BT.bI(end))
                Gs{BT.s+1}=GM*K{indi};
            end
            if logical(BT.bI(end-1))
                %I assume that in all methods considered here, K{BT.s} is different from K{BT.s-1}
                GM=odeIM(t0+dt*BT.cI(indi),K{indi});% Recompute G if necessary
                Gs{BT.s}=GM*K{indi};
            end

        else
            %I assume for each row, there is at least one zero element. in
            %other word, K{indi} must be different from K{indi-1}
            GM=odeIM(t0+dt*BT.cI(indi),K{indi});% Recompute G if necessary
            Gs{indi}=GM*K{indi};

        end
    end

    if BT.is_need_Fs(indi)
        % Check if we need to compute the explicit function for future stages
        Fs{indi}=odeEX(t0+BT.cE(indi)*dt, K{indi});
    end

end


if isfinite(BT.alpha)
    %pt can be computed without adding all K{indi};
    if abs(BT.alpha-1)<20*eps
        pt=K{BT.s};    % Handle special case when alpha is finite and close to 1
    else
        pt=1/BT.alpha*K{BT.s}+(1-1/BT.alpha)*p0;% General case with alpha
    end
else
    % General case where pt is computed as a linear combination of K{indi}

    for indi=1:BT.s
        if logical(BT.bE(indi))
            pt = pt + BT.bE(indi)*Fs{indi};
        end
        if logical(BT.bI(indi))
            pt = pt + BT.bI(indi)*Gs{indi};
        end
    end
    if logical(BT.bI(BT.s+1))
        pt = pt + BT.bI(BT.s+1)*Gs{BT.s+1};
    end

    pt=p0+dt*pt;

end
end