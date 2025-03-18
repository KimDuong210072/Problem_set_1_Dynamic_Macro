%part 3b
function [y, P] = rouwenhorst(N, gamma, sigma, mu)
   sigma_y = sigma / sqrt(1 - gamma^2);
   z = linspace(-sqrt(N-1) * sigma_y, sqrt(N-1) * sigma_y, N);
   y = mu + z; 
   p = (1 + gamma) / 2;
   P = [p, 1 - p; 1 - p, p];
   for n = 3:N
       P_old = P;
       P = p * [P_old, zeros(n-1, 1); zeros(1, n)];
       P = P + (1 - p) * [zeros(n-1, 1), P_old; zeros(1, n)];
       P = P + (1 - p) * [zeros(1, n); P_old, zeros(n-1, 1)];
       P = P + p * [zeros(1, n); zeros(n-1, 1), P_old];
       P(2:end-1, :) = P(2:end-1, :) / 2;
   end
   disp('State Space (y):');
   disp(y);
   disp('Transition Matrix (P):');
   disp(P);
end
N = 7;            
gamma = 0.85;      
sigma = 1;        
mu = 0.5 / (1 - gamma); 
[y, P] = rouwenhorst(N, gamma, sigma, mu)

%%
%part 3c
function [y, P, sim_states] = rouwenhorst_simulation(N, gamma, sigma, mu, T)
   rng(2025);  % Seed 
   sigma_y = sigma / sqrt(1 - gamma^2);
   z = linspace(-sqrt(N-1) * sigma_y, sqrt(N-1) * sigma_y, N);
   y = mu + z;  
   p = (1 + gamma) / 2;
   P = [p, 1 - p; 1 - p, p];
   for n = 3:N
       P_old = P;
       P = p * [P_old, zeros(n-1, 1); zeros(1, n)];
       P = P + (1 - p) * [zeros(n-1, 1), P_old; zeros(1, n)];
       P = P + (1 - p) * [zeros(1, n); P_old, zeros(n-1, 1)];
       P = P + p * [zeros(1, n); zeros(n-1, 1), P_old];
       P(2:end-1, :) = P(2:end-1, :) / 2;
   end
   sim_states = zeros(T, 1);
   sim_states(1) = randi(N); 
      for t = 2:T
       current_state = sim_states(t-1);
       transition_probs = P(current_state, :);
       sim_states(t) = find(mnrnd(1, transition_probs) == 1);
   end
   sim_states = y(sim_states);
   disp('State Space (y):');
   disp(y);
   disp('Transition Matrix (P):');
   disp(P);
   disp('Simulated Markov Chain States:');
   disp(sim_states);
end

rng(2025);  % Set seed for simulation
N = 7;      % Number of states
gamma = 0.85;  
sigma = 1;     
mu = 0.5 / (1 - gamma); 
T = 50; 
[y, P, sim_states] = rouwenhorst_simulation(N, gamma, sigma, mu, T)
figure;
plot(1:T, sim_states, 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('State');
title('Simulated AR(1) Process (Rouwenhorst)');
grid on;
%% Part 3d
rng(2025); 
N = 7;      
gamma = 0.85;  
gamma1 = 0.9;
sigma = 1;     
mu = 0.5 / (1 - gamma); 
T = 50; 
[y, P, sim_states] = rouwenhorst_simulation(N, gamma, sigma, mu, T)
[y, P, sim_states1] = rouwenhorst_simulation(N, gamma1, sigma, mu, T)
tiledlayout(2,2);
plot(1:T, sim_states, 'b', 'LineWidth', 1.5);
title('gamma 1');
hold on;
plot(1:T, sim_states1, 'r', 'LineWidth', 1.5);
title('gamma 2');
hold off;
xlabel('x');
ylabel('y');
title('Gamma');
%title('gamma 2');
%plot(x, y1, 'r', 'LineWidth', 2); % Red sine wave
%hold on;
%plot(x, y2, 'b--', 'LineWidth', 2); 
%figure;
%plot(1:T, sim_states1, 'b', 'LineWidth', 1.5);
%xlabel('Time');
%ylabel('State');
%title('Simulated AR(1) Process (Rouwenhorst)');
grid on;