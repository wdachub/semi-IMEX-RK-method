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