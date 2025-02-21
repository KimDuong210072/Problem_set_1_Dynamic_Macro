%% File Info.

%{

    simulate.m
    ----------
    This code simulates the model.

%}

%% Simulate class.

classdef simulate
    methods(Static)
        %% Simulate the model. 
        
        function sim = grow(par,sol)            
            %% Set up.
            
            kgrid = par.kgrid; % Capital today (state variable).

            yout = sol.y; % Production function.
            kpol = sol.k; % Policy function for capital.
            cpol = sol.c; % Policy function for consumption.
            ipol = sol.i; % Policy function for investment.

            T = par.T; % Time periods.
            ysim = zeros(par.T,1); % Container for simulated output.
            ksim = zeros(par.T,1); % Container for simulated capital stock.
            csim = zeros(par.T,1); % Container for simulated consumption.
            isim = zeros(par.T,1); % Container for simulated investment.
            usim = zeros(par.T,1); % Container for simulated utility.
            
            %% Begin simulation.
            
            rng(par.seed);
            k0_ind = randsample(par.klen,1); % Index for initial capital stock.
            ysim(1) = yout(k0_ind); % Output in period 1 given k0.
            csim(1) = cpol(k0_ind); % Consumption in period 1 given k0.
            ksim(1) = kpol(k0_ind); % Capital choice for period 2 given k0.
            isim(1) = ipol(k0_ind); % Investment in period 1 given k0.
            usim(1) = model.utility(csim(1),par); % Utility in period 1 given k0.

            %% Simulate endogenous variables.

            for j = 2:T % Time loop.
                kt_ind = find(ksim(j-1)==kgrid); % Capital choice in the previous period is the state today. Find where the latter is on the grid.
                ysim(j) = yout(kt_ind); % Output in period t.
                csim(j) = cpol(kt_ind); % Consumption in period t.
                ksim(j) = kpol(kt_ind); % Capital stock for period t+1.
                isim(j) = ipol(kt_ind); % Investment in period t.
                usim(j) = model.utility(csim(j),par); % Utility in period t.
            end

            sim = struct();
            
            sim.ysim = ysim; % Simulated output.
            sim.ksim = ksim; % Simulated capital choice.
            sim.csim = csim; % Simulated consumption.
            sim.isim = isim; % Simulated investment.
            sim.usim = usim; % Simulated utility.
             
        end
        
    end
end