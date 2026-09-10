clc; clear; close all;

%% ------------------------------
%  True First-order models
% ------------------------------
Kp1_true  = 0.1250;     % Model 1 gain
tau1_true = 0.9375;
G1_true = tf(Kp1_true, [tau1_true 1]);

Kp2_true  = 0.5625;       % Model 2 gain
tau2_true = 0.9375;
G2_true = tf(Kp2_true, [tau2_true 1]);

% Time vector
Ts = 1;
t = 0:Ts:300;

%% ------------------------------
% Piecewise input (20% up, 20% down)
% ------------------------------
u = ones(size(t));
u(t > 100) = 1.2;   % +20%
u(t > 200) = 0.8;   % -20%

%% ------------------------------
% System 1 responses
% ------------------------------
y1 = lsim(G1_true, u, t);              % clean output

% Add 2% noise
noise1 = 0.02 * y1 .* randn(size(y1));
y1_noise = y1 + noise1;

%% ------------------------------
% System 2 responses
% ------------------------------
y2 = lsim(G2_true, u, t);              % clean output

% Add 2% noise
noise2 = 0.02 * y2 .* randn(size(y2));
y2_noise = y2 + noise2;

%% ------------------------------
% Plotting real & noisy data
% ------------------------------
figure(1);

% --- System 1 clean ---
subplot(2,2,1);
plot(t, y1, 'LineWidth', 2);    % BLUE
xlabel('Time (s)');
ylabel('y_1(t)');
title('System 1 Clean Output (Kp = 0.133)');
grid on; xlim([0 300]);

% --- System 1 noisy ---
subplot(2,2,2);
plot(t, y1_noise, 'LineWidth', 2);   % RED
xlabel('Time (s)');
ylabel('y_1(t) noisy');
title('System 1 Noisy Output (2% Noise)');
grid on; xlim([0 300]);

% --- System 2 clean ---
subplot(2,2,3);
plot(t, y2, 'LineWidth', 2);    % GREEN
xlabel('Time (s)');
ylabel('y_2(t)');
title('System 2 Clean Output (Kp = 0.6)');
grid on; xlim([0 300]);

% --- System 2 noisy ---
subplot(2,2,4);
plot(t, y2_noise, 'LineWidth', 2);   % MAGENTA
xlabel('Time (s)');
ylabel('y_2(t) noisy');
title('System 2 Noisy Output (2% Noise)');
grid on; xlim([0 300]);


%% ====================================================
% PARAMETER ESTIMATION (WITH LOWER & UPPER BOUNDS)
% Using nonlinear regression with fmincon
% ====================================================

% Cost function
cost_fun = @(x, t, u, y_meas) ...
    sum( (y_meas - lsim(tf(x(1), [x(2) 1]), u, t)).^2 );

% Initial guesses
x0_sys1 = [0.1 1];
x0_sys2 = [0.5 1];

% Lower and upper bounds
lb = [0     0.1];      % [Kp  tau]
ub = [2.0   10];

% fmincon options
opts = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','interior-point', ...
    'MaxIterations',1e4);

%% ----- System 1 (clean data) -----
x_hat1_clean = fmincon(@(x) cost_fun(x, t, u, y1), ...
                       x0_sys1, [], [], [], [], lb, ub, [], opts);
Kp1_est_clean  = x_hat1_clean(1);
tau1_est_clean = x_hat1_clean(2);

%% ----- System 1 (noisy data) -----
x_hat1_noise = fmincon(@(x) cost_fun(x, t, u, y1_noise), ...
                       x0_sys1, [], [], [], [], lb, ub, [], opts);
Kp1_est_noise  = x_hat1_noise(1);
tau1_est_noise = x_hat1_noise(2);

%% ----- System 2 (clean data) -----
x_hat2_clean = fmincon(@(x) cost_fun(x, t, u, y2), ...
                       x0_sys2, [], [], [], [], lb, ub, [], opts);
Kp2_est_clean  = x_hat2_clean(1);
tau2_est_clean = x_hat2_clean(2);

%% ----- System 2 (noisy data) -----
x_hat2_noise = fmincon(@(x) cost_fun(x, t, u, y2_noise), ...
                       x0_sys2, [], [], [], [], lb, ub, [], opts);
Kp2_est_noise  = x_hat2_noise(1);
tau2_est_noise = x_hat2_noise(2);

%% ------------------------------
% Display results
% ------------------------------
fprintf('\n===== TRUE PARAMETERS =====\n');
fprintf('System 1: Kp = %.4f, tau = %.4f\n', Kp1_true, tau1_true);
fprintf('System 2: Kp = %.4f, tau = %.4f\n', Kp2_true, tau2_true);

fprintf('\n===== ESTIMATED FROM CLEAN DATA =====\n');
fprintf('System 1 (clean): Kp = %.4f, tau = %.4f\n', Kp1_est_clean, tau1_est_clean);
fprintf('System 2 (clean): Kp = %.4f, tau = %.4f\n', Kp2_est_clean, tau2_est_clean);

fprintf('\n===== ESTIMATED FROM NOISY DATA =====\n');
fprintf('System 1 (noisy): Kp = %.4f, tau = %.4f\n', Kp1_est_noise, tau1_est_noise);
fprintf('System 2 (noisy): Kp = %.4f, tau = %.4f\n', Kp2_est_noise, tau2_est_noise);



% ------------------------------
% First-order models
% ------------------------------
K1  = 0.1250;     % Model 1 gain
tau1 = 0.9375;
G1 = tf(K1, [tau1 1]);

K2  = 0.5625;       % Model 2 gain
tau2 = 0.9375;
G2 = tf(K2, [tau2 1]);

% Time
Ts = 1;
t = 0:Ts:300;

% Input (+20%, -20%)
u = ones(size(t));
u(t > 100) = 1.2;    % +20%
u(t > 200) = 0.8;    % -20%

% ------------------------------
% System 1: y1 clean & noisy
% ------------------------------
y1 = lsim(G1, u, t);
noise1 = 0.02 * y1 .* randn(size(y1));
y1_noise = y1 + noise1;

% ------------------------------
% System 2: y2 clean & noisy
% ------------------------------
y2 = lsim(G2, u, t);
noise2 = 0.02 * y2 .* randn(size(y2));
y2_noise = y2 + noise2;

% ------------------------------
% PLOTTING (2 GRAPHS)
% ------------------------------
figure(2);

% ===== GRAPH 1: H1 clean vs noisy =====
subplot(2,1,1);
plot(t, y1, 'b', 'LineWidth', 2); hold on;                 % Clean (solid)
plot(t, y1_noise, 'r--', 'LineWidth', 1.8);               % Noisy (dashed)
xlabel('Time (s)'); ylabel('h_1(t)');
title('System H1 : Clean vs Noisy (2% Noise)');
legend('Clean','Noisy (2%)');
grid on; xlim([0 300]);

% ===== GRAPH 2: H2 clean vs noisy =====
subplot(2,1,2);
plot(t, y2, 'b', 'LineWidth', 2); hold on;                 % Clean (solid)
plot(t, y2_noise, 'r-', 'LineWidth', 1.8);               % Noisy (dashed)
xlabel('Time (s)'); ylabel('h_2(t)');
title('System H2 : Clean vs Noisy (2% Noise)');
legend('Clean','Noisy (2%)');
grid on; xlim([0 300]);


