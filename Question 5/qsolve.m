%% Solve class.
classdef qsolve
    methods(Static)
        %% Solve the model using VFI. 
        function sol = grow(par)
            %% Structure array for model solution.
            sol = struct();

            %% Model parameters, grids, and functions.
            beta = par.beta; 
            alpha = par.alpha; 
            gamma = par.gamma; 
            sigma = par.sigma; 
            deltak = par.deltak; 
            deltah = par.deltah; 
            u = par.u; 

            klen = par.klen; 
            kgrid = par.kgrid; 
            hlen = par.hlen; 
            hgrid = par.hgrid; 

            Alen = par.Alen; 
            Agrid = par.Agrid; 
            pmat = par.pmat; 

            % Create 3D grids
            [K, H, A] = ndgrid(kgrid, hgrid, Agrid);
            
            %% Initial guesses for value function and policy functions.
            h0 = H ./ (2 - deltah - u);
            y0 = A .* K.^alpha .* (u .* h0).^gamma; % Ensure y0 is 3D
            i0 = deltak * K; % Ensure i0 is 3D
            c0 = y0 - i0; % Ensure c0 is 3D
            v0 = qmodel.utility(c0, par) ./ (1 - beta); % Ensure v0 is 3D

            v1 = nan(klen, hlen, Alen);
            k1 = nan(klen, hlen, Alen);
            h1 = nan(klen, hlen, Alen);

            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;

            fprintf('\n------------Beginning Value Function Iteration.------------\n\n');

            %% Value Function Iteration
            while diff > crit && iter < maxiter
                for p = 1:klen
                    for s = 1:hlen
                        for j = 1:Alen
                            % Compute output
                            y = Agrid(j) * kgrid(p)^alpha * (u * hgrid(s))^gamma;

                            % Compute investment and consumption (ensure 1D)
                            i = kgrid - (1 - deltak) * kgrid(p);
                            c = y - i;

                            % Compute expected value (Ensure it's a vector)
                            ev = squeeze(v0(:, s, :)) * pmat(j, :)';

                            % Ensure `vall` is a 1D vector
                            vall = qmodel.utility(c, par) + beta * ev;
                            vall = squeeze(vall); % Convert to 1D
                            vall(c <= 0) = -inf; % Handle infeasibility

                            % Get max value & index
                            [vmax, ind] = max(vall(:));

                            % Assign scalar values correctly
                            v1(p, s, j) = vmax;
                            k1(p, s, j) = kgrid(ind);
                            h1(p, s, j) = hgrid(ind);
                        end
                    end
                end

                % Compute norm using Frobenius norm for 3D matrices
                diff = norm(v1 - v0, 'fro');  
                v0 = v1;
                iter = iter + 1;

                % Mark every 25 iterations for visibility
                if mod(iter,25) == 0
                    fprintf('Iteration: %d.\n', iter);
                end
            end

            fprintf('\nConverged in %d iterations.\n', iter);
            fprintf('\n------------End of Value Function Iteration.------------\n\n');

            %% Macro variables and policy functions
            sol.y = y0;  
            sol.k = k1;
            sol.h = h1;
            sol.i = k1 - (1 - deltak) * K; % Ensure i is 3D
            sol.c = sol.y - sol.i;
            sol.v = v1;
        end
    end
end
