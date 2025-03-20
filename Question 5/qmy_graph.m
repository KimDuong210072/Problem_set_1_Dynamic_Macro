%% File Info.
%{
    my_graph.m
    ----------
    This code plots the value and policy functions in 3D.
%}

%% Graph class.
classdef qmy_graph
    methods(Static)
        %% Plot value and policy functions.
        function [] = plot_policy(par,sol,sim,figout)
            %% Select an index for Agrid (Productivity level)
            Aidx = round(par.Alen / 2); % Choose the middle index for visualization

            %% Plot production function.
            figure(1)
            surf(par.kgrid, par.hgrid, sol.y(:,:,Aidx)')
            xlabel({'$k_{t}$'},'Interpreter','latex')
            ylabel({'$h_{t}$'},'Interpreter','latex')
            zlabel({'$y_{t}$'},'Interpreter','latex')
            title('Production Function')
            
            %fig_name = strcat(figout, 'ypol.fig');
            %savefig(fig_name)

            %% Plot capital policy function.
            figure(2)
            surf(par.kgrid, par.hgrid, sol.k(:,:,Aidx)')
            xlabel({'$k_{t}$'},'Interpreter','latex')
            ylabel({'$h_{t}$'},'Interpreter','latex')
            zlabel({'$k_{t+1}$'},'Interpreter','latex')
            title('Capital Policy Function')
            
            %fig_name = strcat(figout, 'kpol.fig');
            %savefig(fig_name)

            %% Plot consumption policy function.
            figure(3)
            surf(par.kgrid, par.hgrid, sol.c(:,:,Aidx)')
            xlabel({'$k_{t}$'},'Interpreter','latex')
            ylabel({'$h_{t}$'},'Interpreter','latex')
            zlabel({'$c_{t}$'},'Interpreter','latex')
            title('Consumption Policy Function')
            
            %fig_name = strcat(figout, 'cpol.fig');
            %savefig(fig_name)

            %% Plot investment policy function.
            figure(4)
            surf(par.kgrid, par.hgrid, sol.i(:,:,Aidx)')
            xlabel({'$k_{t}$'},'Interpreter','latex')
            ylabel({'$h_{t}$'},'Interpreter','latex')
            zlabel({'$i_{t}$'},'Interpreter','latex')
            title('Investment Policy Function')
            
            %fig_name = strcat(figout, 'ipol.fig');
            %savefig(fig_name)

            %% Plot value function.
            figure(5)
            surf(par.kgrid, par.hgrid, sol.v(:,:,Aidx)')
            xlabel({'$k_{t}$'},'Interpreter','latex')
            ylabel({'$h_{t}$'},'Interpreter','latex')
            zlabel({'$v_t(k_t,h_t,A_t)$'},'Interpreter','latex')
            title('Value Function')
            
            %fig_name = strcat(figout, 'vfun.fig');
            %savefig(fig_name)

            %% Plot simulated output.
            tgrid = linspace(1, par.T, par.T);
            figure(6)
            plot(tgrid, sim.ysim)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$y^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Output')
            %fig_name = strcat(figout, 'ysim.fig');
            %savefig(fig_name)

            %% Plot simulated capital choice.
            figure(7)
            plot(tgrid, sim.ksim)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$k^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Capital Choice')
            %fig_name = strcat(figout, 'ksim.fig');
            %savefig(fig_name)

            %% Plot simulated consumption.
            figure(8)
            plot(tgrid, sim.csim)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$c^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Consumption')
            %fig_name = strcat(figout, 'csim.fig');
            %savefig(fig_name)

            %% Plot simulated investment.
            figure(9)
            plot(tgrid, sim.isim)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$i^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Investment')
            %fig_name = strcat(figout, 'isim.fig');
            %savefig(fig_name)

            %% Plot simulated utility.
            figure(10)
            plot(tgrid, sim.usim)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$u^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Utility')
            %fig_name = strcat(figout, 'usim.fig');
            %savefig(fig_name)

            %% Plot simulated productivity.
            figure(11)
            plot(tgrid, sim.Asim)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$A^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Productivity')
            %fig_name = strcat(figout, 'Asim.fig');
            %savefig(fig_name)
        end
    end
end
