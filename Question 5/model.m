%% File Info.

%{

    model.m
    -------
    This code sets up the model.

%}

%% Model class.

classdef model
    methods(Static)
        %% Set up structure array for model parameters and set the simulation parameters.
        
        function par = setup()            
            %% Structure array for model parameters.
            
            par = struct();
            
            %% Preferences.
            
            par.beta = 0.96; % Discount factor: Lower values of this mean that consumers are impatient and consume more today.
            par.sigma = 7.0; % CRRA: Higher values of this mean that consumers are risk averse and do not want to consume too much today.
            
            assert(par.beta > 0.0 && par.beta < 1.0,'Discount factor should be between 0 and 1.\n')
            assert(par.sigma > 0.0,'CRRA should be at least 0.\n')

            %% Technology.

            par.alpha = 0.4; % Physical capital's share of income.
            par.gamma = 1-par.alpha; % Human capital's share of income.
            par.delta_k = 0.08; % Depreciation rate of physical capital.
            par.delta_h = 0.6; % Depreciation rate of human capital.

            assert(par.alpha > 0.0 && par.alpha < 1.0,'Physical capital share of income should be between 0 and 1.\n')
            assert(par.gamma > 0.0 && par.gamma < 1.0,'Human capital share of income should be between 0 and 1.\n')
            assert(par.alpha + par.gamma == 1.0 ,'Production should have constant returns to scale.\n')
            assert(par.delta_k >= 0.0 && par.delta_k <= 1.0,'The depreciation rate should be from 0 to 1.\n')
            assert(par.delta_h >= 0.0 && par.delta_h <= 1.0,'The depreciation rate should be from 0 to 1.\n')

            par.sigma_eps = 0.09; % Std. dev of productivity shocks.
            par.rho = 0.87; % Persistence of AR(1) process.
            par.mu = 0.0; % Intercept of AR(1) process.

            assert(par.sigma_eps > 0.0,'The standard deviation of the shock must be positive.\n')
            assert(abs(par.rho) < 1.0,'The persistence must be less than 1 in absolute value so that the series is stationary.\n')

            %% Simulation parameters.

            par.seed = 2025; % Seed for simulation.
            par.T = 20; % Number of time periods.

        end
        
        %% Generate state grids.
        
        function par = gen_grids(par)
            %% Capital grid.

            par.klen = 20; % Grid size for k.
            par.kmax = 10.0; % Upper bound for k.
            par.kmin = 2.0; % Minimum k.
            
            assert(par.klen > 5,'Grid size for k should be positive and greater than 5.\n')
            assert(par.kmax > par.kmin,'Minimum k should be less than maximum value.\n')
            assert(par.kmin > 0.0,'Capital should always be positive.\n')
            
            par.kgrid = linspace(par.kmin,par.kmax,par.klen)'; % Equally spaced, linear grid for k and k'.
            
            %% Human capital grid.

            par.hlen = 20; % Grid size for h.
            par.hmax = 10.0; % Upper bound for h.
            par.hmin = 2.0; % Minimum h.
            
            assert(par.hlen > 5,'Grid size for h should be positive and greater than 5.\n')
            assert(par.hmax > par.hmin,'Minimum h should be less than maximum value.\n')
            assert(par.hmin > 0.0,'Human capital should always be positive.\n')
            
            par.hgrid = linspace(par.hmin,par.hmax,par.hlen)'; % Equally spaced, linear grid for h and h'.

            %% Discretized productivity process.
                  
            par.Alen = 10; % Grid size for A.
            par.m = 3; % Scaling parameter for Tauchen.
            
            assert(par.Alen > 3,'Grid size for A should be positive and greater than 3.\n')
            assert(par.m > 0,'Scaling parameter for Tauchen should be positive.\n')
            
            [Agrid,pmat] = model.tauchen(par.mu,par.rho,par.sigma_eps,par.Alen,par.m); % Tauchen's Method to discretize the AR(1) process for log productivity.
            par.Agrid = exp(Agrid); % The AR(1) is in logs so exponentiate it to get A.
            par.pmat = pmat; % Transition matrix.

        end
        
        %% Tauchen's Method
        
        function [y,pi] = tauchen(mu,rho,sigma,N,m)
            %% Construct equally spaced grid.
        
            ar_mean = mu/(1-rho); % The mean of a stationary AR(1) process is mu/(1-rho).
            ar_sd = sigma/((1-rho^2)^(1/2)); % The std. dev of a stationary AR(1) process is sigma/sqrt(1-rho^2)
            
            y1 = ar_mean-(m*ar_sd); % Smallest grid point is the mean of the AR(1) process minus m*std.dev of AR(1) process.
            yn = ar_mean+(m*ar_sd); % Largest grid point is the mean of the AR(1) process plus m*std.dev of AR(1) process.
            
	        y = linspace(y1,yn,N); % Equally spaced grid.
            d = y(2)-y(1); % Step size.
	        
	        %% Compute transition probability matrix from state j (row) to k (column).
        
            ymatk = repmat(y,N,1); % States next period.
            ymatj = mu+rho*ymatk'; % States this period.
        
	        pi = normcdf(ymatk,ymatj-(d/2),sigma) - normcdf(ymatk,ymatj+(d/2),sigma); % Transition probabilities to state 2, ..., N-1.
	        pi(:,1) = normcdf(y(1),mu+rho*y-(d/2),sigma); % Transition probabilities to state 1.
	        pi(:,N) = 1 - normcdf(y(N),mu+rho*y+(d/2),sigma); % Transition probabilities to state N.
	        
        end
        
        %% Utility function.
        
        function u = utility(c,par)
            %% CRRA utility.
            
            if par.sigma == 1
                u = log(c); % Log utility.
            else
                u = (c.^(1-par.sigma))./(1-par.sigma); % CRRA utility.
            end
                        
        end
        
    end
end