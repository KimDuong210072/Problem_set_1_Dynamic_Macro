%% File Info.

%{

    solve.m
    -------
    This code solves the model.

%}

%% Solve class.

classdef solve
    methods(Static)
        %% Solve the model using VFI. 
        
        function sol = grow(par)            
            %% Structure array for model solution.
            
            sol = struct();
            
            %% Model parameters, grids and functions.
            
            beta = par.beta; % Discount factor.
            alpha = par.alpha; % Capital share of income.
            delta = par.delta; % Depreciation rate.

            klen = par.klen; % Grid size for k.
            kgrid = par.kgrid; % Grid for k (state and choice).

            %% Value Function Iteration.
            
            y0 = kgrid.^alpha; % Output in the "last period", given k.
            i0 = (1-delta).*kgrid; % In the "last period," k' is zero but k still depreciates.
            c0 = y0-i0; % Consumption in "last period."
            v0 = model.utility(c0,par)./(1-beta); % Guess of value function for each value of k.

            v1 = zeros(klen,1); % Container for V.
            k1 = zeros(klen,1); % Container for k'.
                            
            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;
            
            fprintf('------------Beginning Value Function Iteration.------------\n\n')
            
            while diff > crit && iter < maxiter % Iterate on the Bellman Equation until convergence.
                
                for p = 1:klen % Loop over the k-states.
                    
                    % Macro variables.
                    y = kgrid(p)^alpha; % Output given k, kgrid(p).
                    i = kgrid-((1-delta)*kgrid(p)); % Possible values for investment, i=k'-(1-delta)k, when choosing k' from kgrid and given k.
                    c = y-i; % Possible values for consumption, c = y-i, given y and i.
                    
                    % Solve the maximization problem.
                    vall = model.utility(c,par) + beta*v0; % Compute the value function for each choice of k', given k.
                    vall(c<0) = -inf; % Set the value function to negative infinity when c < 0.
                    [vmax,ind] = max(vall); % Maximize: vmax is the maximized value function; ind is where it is in the grid.
                    
                    % Store values.
                    v1(p) = vmax; % Maximized v.
                    k1(p) = kgrid(ind); % Optimal k'.

                end
                
                diff = norm(v1-v0); % Check for convergence.
                v0 = v1; % Update guess of v.
                
                iter = iter + 1; % Update counter.
                
                % Print counter.
                if mod(iter,25) == 0
                    fprintf('Iteration: %d.\n',iter)
                end

            end
                
            fprintf('\nConverged in %d iterations.\n\n',iter)
            
            fprintf('------------End of Value Function Iteration.------------\n')
            
            %% Macro variables, value, and policy functions.
            
            sol.y = kgrid.^alpha; % Output.
            sol.k = k1; % Capital policy function.
            sol.i = k1-((1-delta).*kgrid); % Investment policy function.
            sol.c = sol.y-sol.i; % Consumption policy function.
            sol.v = v1; % Value function.
            
        end
        
    end
end