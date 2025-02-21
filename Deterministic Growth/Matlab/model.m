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
            par.sigma = 2.00; % CRRA: Higher values of this mean that consumers are risk averse and do not want to consume too much today.
            
            assert(par.beta > 0 && par.beta < 1.00,'Discount factor should be between 0 and 1.\n')
            assert(par.sigma > 0,'CRRA should be at least 0.\n')

            %% Technology.

            par.alpha = 0.33; % Capital's share of income.
            par.delta = 0.05; % Depreciation rate of physical capital.

            assert(par.alpha > 0 && par.alpha < 1.00,'Capital share of income should be between 0 and 1.\n')
            assert(par.delta >= 0 && par.delta <= 1.00,'The depreciation rate should be from 0 to 1.\n')

            %% Simulation parameters.

            par.seed = 2025; % Seed for simulation.
            par.T = 100; % Number of time periods.

        end
        
        %% Generate state grids.
        
        function par = gen_grids(par)
            %% Capital grid.

            par.kss = (par.alpha/((1/par.beta)-1+par.delta))^(1/(1-par.alpha)); % Steady state capital.
             
            par.klen = 300; % Grid size for k.
            par.kmax = 1.75*par.kss; % Upper bound for k.
            par.kmin = 0.25*par.kss; % Minimum k.
            
            assert(par.klen > 5,'Grid size for k should be positive and greater than 5.\n')
            assert(par.kmax > par.kmin,'Minimum k should be less than maximum value.\n')
            
            par.kgrid = linspace(par.kmin,par.kmax,par.klen)'; % Equally spaced, linear grid for k and k'.
                        
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