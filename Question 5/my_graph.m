%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions.

%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_policy(par,sol,sim)
            %% Plot production function.
            
            figure(1)

            surf(par.kgrid,par.hgrid,squeeze(sol.y(:,1,:))')
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$h_{t}$'},'Interpreter','latex') 
                zlabel({'$y_{t}$'},'Interpreter','latex') 
            title('Production Function')

            %% Plot capital policy function.
            
            figure(2)
            
            surf(par.kgrid,par.hgrid,squeeze(sol.k(:,1,:))')
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$h_{t}$'},'Interpreter','latex') 
                zlabel({'$k_{t+1}$'},'Interpreter','latex') 
            title('Capital Policy Function')
            
            %% Plot human capital policy function.
            
            figure(3)
            
            surf(par.kgrid,par.hgrid,squeeze(sol.h(:,1,:))')
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$h_{t}$'},'Interpreter','latex') 
                zlabel({'$h_{t+1}$'},'Interpreter','latex') 
            title('Human Capital Policy Function')
            
            %% Plot consumption policy function.
            
            figure(4)
            
            surf(par.kgrid,par.hgrid,squeeze(sol.c(:,1,:))')
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$h_{t}$'},'Interpreter','latex') 
                zlabel({'$c_{t}$'},'Interpreter','latex') 
            title('Consumption Policy Function')
            
            %% Plot capital investment policy function.
            
            figure(5)
            
            surf(par.kgrid,par.hgrid,squeeze(sol.i(:,1,:))')
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$h_{t}$'},'Interpreter','latex') 
                zlabel({'$i^k_{t}$'},'Interpreter','latex') 
            title('Capital Investment Policy Function')
            
            %% Plot human capital investment policy function.
            
            figure(6)
            
            surf(par.kgrid,par.hgrid,squeeze(sol.u(:,1,:))')
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$h_{t}$'},'Interpreter','latex') 
                zlabel({'$1-u_{t}$'},'Interpreter','latex') 
            title('Human Capital Investment Policy Function')
            
            %% Plot value function.
            
            figure(7)
            
            surf(par.kgrid,par.hgrid,squeeze(sol.v(:,1,:))')
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$h_{t}$'},'Interpreter','latex') 
                zlabel({'$v_t(k_t,A_t)$'},'Interpreter','latex') 
            title('Value Function')
            
            %% Plot simulated output.

            tgrid = linspace(1,par.T,par.T);

            figure(8)

            plot(tgrid,sim.ysim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$y^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Output')

            %% Plot simulated capital choice.

            figure(9)

            plot(tgrid,sim.ksim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$k^{sim}_{t+1}$'},'Interpreter','latex') 
            title('Simulated Capital Choice')

            %% Plot simulated capital choice.

            figure(10)

            plot(tgrid,sim.hsim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$h^{sim}_{t+1}$'},'Interpreter','latex') 
            title('Simulated Human Capital Choice')

            %% Plot simulated consumption.

            figure(11)

            plot(tgrid,sim.csim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$c^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Consumption')

            %% Plot simulated capital investment.

            figure(12)

            plot(tgrid,sim.isim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$i^{ksim}_t$'},'Interpreter','latex') 
            title('Simulated Capital Investment')

            %% Plot simulated human capital investment.

            figure(13)

            plot(tgrid,sim.usim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$1-u_t$'},'Interpreter','latex') 
            title('Simulated Human Capital Investment')

            %% Plot simulated utility.

            figure(14)

            plot(tgrid,sim.vsim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$v^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Utility')

            %% Plot simulated productivity.

            figure(15)

            plot(tgrid,sim.Asim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$A^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Productivity')

        end
        
    end
end