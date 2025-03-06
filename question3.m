%part 3b
function [y, P] = rouwenhorst(N, gamma, sigma, mu)
    % Rouwenhorst's method to approximate an AR(1) process
    % Inputs:
    %   N     - Number of states
    %   gamma - AR(1) coefficient
    %   sigma - Standard deviation of shocks
    %   mu    - Mean of the process
    % Outputs:
    %   y     - State space (discretized values)
    %   P     - Transition probability matrix

    % Compute standard deviation of stationary distribution
    sigma_y = sigma / sqrt(1 - gamma^2);

    % Define state space (equally spaced points)
    z = linspace(-sqrt(N-1) * sigma_y, sqrt(N-1) * sigma_y, N);
    y = mu + z;  % Shift by mean

    % Compute probability p
    p = (1 + gamma) / 2;

    % Base case for N = 2
    P = [p, 1 - p; 1 - p, p];

    % Recursively construct for N > 2
    for n = 3:N
        P_old = P;
        P = p * [P_old, zeros(n-1, 1); zeros(1, n)];
        P = P + (1 - p) * [zeros(n-1, 1), P_old; zeros(1, n)];
        P = P + (1 - p) * [zeros(1, n); P_old, zeros(n-1, 1)];
        P = P + p * [zeros(1, n); zeros(n-1, 1), P_old];

        % Normalize middle rows
        P(2:end-1, :) = P(2:end-1, :) / 2;
    end

    % Display results
    disp('State Space (y):');
    disp(y);
    disp('Transition Matrix (P):');
    disp(P);
end

% Define parameters
N = 7;             % Number of states
gamma = 0.85;      % Persistence parameter
sigma = 1;         % Std. deviation of shocks
mu = 0.5 / (1 - gamma); % Long-run mean

% Call function
[y, P] = rouwenhorst(N, gamma, sigma, mu);
%%
%part 3c
function [y, P, sim_states] = rouwenhorst_simulation(N, gamma, sigma, mu, T)
    % Set seed for reproducibility
    rng(2025);  % Seed for simulation

    % Compute standard deviation of stationary distribution
    sigma_y = sigma / sqrt(1 - gamma^2);

    % Define state space (evenly spaced grid)
    z = linspace(-sqrt(N-1) * sigma_y, sqrt(N-1) * sigma_y, N);
    y = mu + z;  % Shift by mean

    % Compute probability p
    p = (1 + gamma) / 2;

    % Base case for N = 2
    P = [p, 1 - p; 1 - p, p];

    % Recursively construct for N > 2
    for n = 3:N
        P_old = P;
        P = p * [P_old, zeros(n-1, 1); zeros(1, n)];
        P = P + (1 - p) * [zeros(n-1, 1), P_old; zeros(1, n)];
        P = P + (1 - p) * [zeros(1, n); P_old, zeros(n-1, 1)];
        P = P + p * [zeros(1, n); zeros(n-1, 1), P_old];

        % Normalize middle rows
        P(2:end-1, :) = P(2:end-1, :) / 2;
    end

    % Simulate Markov chain
    sim_states = zeros(T, 1);

    % Step 1: Draw the initial state index assuming a uniform distribution
    sim_states(1) = randi(N); % Random index from 1 to N (Uniform)

    % Step 2: Transition through time using probabilities
    for t = 2:T
        current_state = sim_states(t-1);
        transition_probs = P(current_state, :);
        sim_states(t) = find(mnrnd(1, transition_probs) == 1); % Sample next state
    end

    % Convert state indices to state values
    sim_states = y(sim_states);

    % Display results
    disp('State Space (y):');
    disp(y);
    disp('Transition Matrix (P):');
    disp(P);
    disp('Simulated Markov Chain States:');
    disp(sim_states);
end

% Define parameters
rng(2025);  % Set seed for simulation
N = 100;      % Number of states
gamma = 0.85;  % Persistence parameter
sigma = 1;     % Std. deviation of shocks
mu = 0.5 / (1 - gamma); % Long-run mean
T = 100;    % Number of periods to simulate

% Call function
[y, P, sim_states] = rouwenhorst_simulation(N, gamma, sigma, mu, T)

figure;
plot(1:T, sim_states, 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('State');
title('Simulated AR(1) Process (Rouwenhorst)');
grid on;
