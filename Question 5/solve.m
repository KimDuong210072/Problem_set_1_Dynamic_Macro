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
            gamma = par.gamma; % Human capital share of income.
            delta_k = par.delta_k; % Depreciation rate (capital).
            delta_h = par.delta_h; % Depreciation rate (human capital).

            klen = par.klen; % Grid size for k.
            kgrid = par.kgrid; % Grid for k (state and choice).
            hlen = par.hlen; % Grid size for h.
            hgrid = par.hgrid; % Grid for h (state and choice).

            Alen = par.Alen; % Grid size for A.
            Agrid = par.Agrid; % Grid for A.
            pmat = par.pmat; % Grid for A.

            kmat = repmat(kgrid,1,Alen,hlen); % k for each value of A and h.
            Amat = repmat(Agrid,klen,1,hlen); % A for each value of k and h.
            hmat = repmat(reshape(hgrid,1,1,[]),klen,Alen,1); % h for each value of k and A.

            %% Value Function Iteration.
            
            y0 = Amat.*kmat.^alpha.*((2-delta_h).*hmat).^gamma; % Fix the inputs. Assume h'=0.
            i0 = -(1-delta_k)*kmat; % Assume k'=0.
            c0 = y0-i0; % Steady-state consumption in a deterministic setting.
            v0 = model.utility(c0,par)./(1-beta); % Guess of value function for each combination of k and A; each row corresponds to a given k and each column corresponds to a given A.

            v1 = zeros(klen,Alen,hlen); % Container for V.
            k1 = zeros(klen,Alen,hlen); % Container for k'.
            h1 = zeros(klen,Alen,hlen); % Container for h'.
            u1 = zeros(klen,Alen,hlen); % Container for u.
            y1 = zeros(klen,Alen,hlen); % Container for y.
            i1 = zeros(klen,Alen,hlen); % Container for i.
            c1 = zeros(klen,Alen,hlen); % Container for c.
                            
            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;
            
            fprintf('------------Beginning Value Function Iteration.------------\n\n')
            
            while diff > crit && iter < maxiter % Iterate on the Bellman Equation until convergence.
                
                for j = 1:Alen % Loop over the A-states.
                    
                    % Pre-compute the expectation.
                    ev = zeros(klen,hlen);
                    for pp = 1:klen
                        for qq = 1:hlen
                           ev(pp,qq) = squeeze(v0(pp,:,qq))*pmat(j,:)'; % The next-period value function is the expected value function over each possible next-period A, conditional on the current state j.
                        end
                    end

                    for p = 1:klen % Loop over the k-states.
                        for q = 1:hlen % Loop over the h-states.
                            
                            % Choices for u, the fraction of time spent on production.
                            u = 2.0-delta_h-(hgrid'./hgrid(q)); % Possible values for u, u=2-delta_h-(h'/h), when choosing h' from hgrid and given h. Transpose u so c is a matrix whose dimensions are klen by hlen.

                            % Find h' such that u > 0.
                            hub = find(u>0.0); % Find the largest value of h' so that h'/h will give u > 0.
                            hub = hub(end);

                            % Find h' such that u <= 1.
                            hlb = find(u<=1.0); % Find the smallest value of h' so that h'/h will give u <= 1.
                            hlb = hlb(1);

                            % Adjust hgrid for u that are out of bounds.
                            hp = hgrid; 
                            hp(u>1.0) = hgrid(hlb);
                            hp(u<=0.0) = hgrid(hub);
                            u = 2.0-delta_h-(hp'./hgrid(q)); % Possible values for u, u=2-delta_h-(h'/h), when choosing h' from hgrid and given h. Transpose u so c is a matrix whose dimensions are klen by hlen.

                            % Macro variables.
                            y = Agrid(j)*kgrid(p)^alpha*(u.*hgrid(q)).^gamma; % Output given k, h, and A, kgrid(p), hgrid(q), u, and Agrid(j) respectively.
                            i = kgrid-(1-delta_k)*kgrid(p); % Possible values for investment, i_k=k'-(1-delta_k)k, when choosing k' from kgrid and given k.
                            c = y-i; % Possible values for consumption, c = y-i, given y and i.
                            c(c<0.0) = 0.0;
    
                            % Solve the maximization problem.
                            vall = model.utility(c,par) + beta*ev; % Compute the value function for each choice of k' and h', given k and h respectively.
                            vall(c<=0.0) = -inf; % Set the value function to negative infinity when c < 0.
                            [vmax,ind] = max(vall,[],'all','linear'); % Maximize: vmax is the maximized value function; ind is where it is in the grid.
                            [kind,hind] = ind2sub([klen hlen],ind);

                            % Store values.
                            v1(p,j,q) = vmax; % Maximized v.
                            k1(p,j,q) = kgrid(kind); % Optimal k'.
                            u1(p,j,q) = 1-u(hind); % Optimal u.
                            h1(p,j,q) = hp(hind); % Optimal h'.
                            i1(p,j,q) = i(kind); % Optimal i_k.
                            c1(p,j,q) = c(kind,hind); % Optimal c.
                            y1(p,j,q) = y(hind); % Optimal y.
    
                        end
                    end
                end
                
                diff = norm(v1-v0,"fro"); % Check for convergence.
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
            
            sol.y = y1; % Output.
            sol.k = k1; % Capital policy function.
            sol.h = h1; % Human capital policy function.
            sol.i = i1; % Capital investment policy function.
            sol.u = u1; % Human capital investment policy function; this is 1-u_t.
            sol.c = c1; % Consumption policy function.
            sol.v = v1; % Value function.
            
        end
        
    end
end