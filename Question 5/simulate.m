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
            hgrid = par.hgrid; % Human capital today (state variable).
            Agrid = par.Agrid; % Productivity (state variable).

            yout = sol.y; % Production function.
            kpol = sol.k; % Policy function for capital.
            hpol = sol.h; % Policy function for human capital.
            cpol = sol.c; % Policy function for consumption.
            ipol = sol.i; % Policy function for investment.
            upol = sol.u; % Policy function for human capital investment; this is 1-u_t.

            T = par.T; % Time periods.
            Asim = zeros(par.T,1); % Container for simulated productivity.
            ysim = zeros(par.T,1); % Container for simulated output.
            ksim = zeros(par.T,1); % Container for simulated capital stock.
            hsim = zeros(par.T,1); % Container for simulated human capital stock.
            csim = zeros(par.T,1); % Container for simulated consumption.
            isim = zeros(par.T,1); % Container for simulated capital investment.
            usim = zeros(par.T,1); % Container for simulated human capital investment, 1-u_t.
            vsim = zeros(par.T,1); % Container for simulated utility.
            
            %% Begin simulation.
            
            rng(par.seed);

            pmat0 = par.pmat^1000;
            pmat0 = pmat0(1,:); % Stationary distribution.
            cmat = cumsum(par.pmat,2); % CDF matrix.

            k0_ind = 1; %randsample(par.klen,1); % Index for initial capital stock.
            h0_ind = 1; %randsample(par.hlen,1); % Index for initial human capital stock.
            A0_ind = randsample(par.Alen,1,true,pmat0); % Index for initial productivity.

            Asim(1) = Agrid(A0_ind); % Productivity in period 1.
            ysim(1) = yout(k0_ind,A0_ind,h0_ind); % Output in period 1 given k0 and A0.
            csim(1) = cpol(k0_ind,A0_ind,h0_ind); % Consumption in period 1 given k0 and A0.
            ksim(1) = kpol(k0_ind,A0_ind,h0_ind); % Capital choice for period 2 given k0 and A0.
            hsim(1) = hpol(k0_ind,A0_ind,h0_ind); % Capital choice for period 2 given k0 and A0.
            isim(1) = ipol(k0_ind,A0_ind,h0_ind); % Capital investment in period 1 given k0 and A0.
            usim(1) = upol(k0_ind,A0_ind,h0_ind); % Human capital investment, 1-u_t, in period 1 given k0 and A0.
            vsim(1) = model.utility(csim(1),par); % Utility in period 1 given k0 and A0.

            A1_ind = find(rand<=cmat(A0_ind,:)); % Draw productivity for next period.
            A0_ind = A1_ind(1);

            %% Simulate endogenous and exogenous variables.

            for j = 2:T % Time loop.
                kt_ind = find(ksim(j-1)==kgrid); % Capital choice in the previous period is the state today. Find where the latter is on the grid.
                ht_ind = find(hsim(j-1)==hgrid); % Human capital choice in the previous period is the state today. Find where the latter is on the grid.
                Asim(j) = Agrid(A0_ind); % Productivity in period t.
                ysim(j) = yout(kt_ind,A0_ind,ht_ind); % Output in period t.
                csim(j) = cpol(kt_ind,A0_ind,ht_ind); % Consumption in period t.
                hsim(j) = hpol(kt_ind,A0_ind,ht_ind); % Human capital in period t+1.
                ksim(j) = kpol(kt_ind,A0_ind,ht_ind); % Capital stock for period t+1.
                isim(j) = ipol(kt_ind,A0_ind,ht_ind); % Capital investment in period t.
                usim(j) = upol(kt_ind,A0_ind,ht_ind); % Human capital investment, 1-u_t, in period t.
                vsim(j) = model.utility(csim(j),par); % Utility in period t.
                A1_ind = find(rand<=cmat(A0_ind,:)); % Draw next state.
                A0_ind = A1_ind(1); % State next period.
            end

            sim = struct();
            
            % Burn the first half.
            sim.Asim = Asim; % Simulated productivity.
            sim.ysim = ysim; % Simulated output.
            sim.ksim = ksim; % Simulated capital choice.
            sim.hsim = hsim; % Simulated human capital choice.
            sim.csim = csim; % Simulated consumption.
            sim.isim = isim; % Simulated capital investment.
            sim.usim = usim; % Simulated human capital investment, 1-u_t.
            sim.vsim = vsim; % Simulated utility.
             
        end
        
    end
end