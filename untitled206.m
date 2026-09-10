clc; clear; close all;

%% ====================================================
%  TRUE FIRST-ORDER (FOPTD) MODELS
% ====================================================
Kp1_true  = 0.1250;
tau1_true = 0.9375;
G1_true = tf(Kp1_true, [tau1_true 1]);

Kp2_true  = 0.5625;
tau2_true = 0.9375;
G2_true = tf(Kp2_true, [tau2_true 1]);

theta_true = 0;   % negligible dead time

%% ====================================================
%  INPUT SIGNAL
% ====================================================
Ts = 1;
t = 0:Ts:300;

u = ones(size(t));
u(t > 100) = 1.2;   % +20%
u(t > 200) = 0.8;   % -20%

%% ====================================================
%  GENERATE PROCESS DATA (CLEAN + NOISY)
% ====================================================
y1 = lsim(G1_true, u, t);
y2 = lsim(G2_true, u, t);

noise_level = 0.02;     % 2% noise
y1_noise = y1 .* (1 + noise_level*randn(size(y1)));
y2_noise = y2 .* (1 + noise_level*randn(size(y2)));

%% ====================================================
%  RAW DATA PLOTS (CLEAN vs NOISY)
% ====================================================
figure(1);

subplot(2,2,1);
plot(t, y1, 'LineWidth', 2);
title('System 1 Clean Output');
xlabel('Time'); ylabel('h_1');
grid on

subplot(2,2,2);
plot(t, y1_noise, 'LineWidth', 2);
title('System 1 Noisy Output (2%)');
xlabel('Time'); ylabel('h_1');
grid on

subplot(2,2,3);
plot(t, y2, 'LineWidth', 2);
title('System 2 Clean Output');
xlabel('Time'); ylabel('h_2');
grid on

subplot(2,2,4);
plot(t, y2_noise, 'LineWidth', 2);
title('System 2 Noisy Output (2%)');
xlabel('Time'); ylabel('h_2');
grid on

%% ====================================================
%  TRAINING / VALIDATION SPLIT (70% / 30%)
% ====================================================
N = length(t);
N_train = round(0.7 * N);

t_train = t(1:N_train);
u_train = u(1:N_train);
y1_train = y1_noise(1:N_train);
y2_train = y2_noise(1:N_train);

t_val = t(N_train+1:end);
u_val = u(N_train+1:end);
y1_val = y1_noise(N_train+1:end);
y2_val = y2_noise(N_train+1:end);

%% ====================================================
%  PARAMETER ESTIMATION USING NONLINEAR REGRESSION
% ====================================================
cost_fun = @(x, t, u, y) ...
    sum((y - lsim(tf(x(1), [x(2) 1]), u, t)).^2);

opts = optimset('Display','iter','MaxIter',1e4,'MaxFunEvals',1e4);

% ----- System 1 (h1) -----
x0_1 = [0.1 1];
xhat1 = fminsearch(@(x) cost_fun(x, t_train, u_train, y1_train), x0_1, opts);
Kp1_id  = xhat1(1);
tau1_id = xhat1(2);

% ----- System 2 (h2) -----
x0_2 = [0.5 1];
xhat2 = fminsearch(@(x) cost_fun(x, t_train, u_train, y2_train), x0_2, opts);
Kp2_id  = xhat2(1);
tau2_id = xhat2(2);

theta_id = 0;   % dead time negligible

%% ====================================================
%  PARAMETER COMPARISON (PART i)
% ====================================================
fprintf('\n============ PARAMETER COMPARISON ============\n');

fprintf('\nSystem 1 (h1/qi)\n');
fprintf('True      : Kp = %.4f , Tau = %.4f\n',Kp1_true,tau1_true);
fprintf('Identified: Kp = %.4f , Tau = %.4f , Theta = %.4f\n',Kp1_id,tau1_id,theta_id);

fprintf('\nSystem 2 (h2/qi)\n');
fprintf('True      : Kp = %.4f , Tau = %.4f\n',Kp2_true,tau2_true);
fprintf('Identified: Kp = %.4f , Tau = %.4f , Theta = %.4f\n',Kp2_id,tau2_id,theta_id);

%% ====================================================
%  VALIDATION USING 30% UNSEEN DATA (PART ii)
% ====================================================
G1_id = tf(Kp1_id, [tau1_id 1], 'InputDelay', theta_id);
G2_id = tf(Kp2_id, [tau2_id 1], 'InputDelay', theta_id);

y1_val_pred = lsim(G1_id, u_val, t_val);
y2_val_pred = lsim(G2_id, u_val, t_val);

%% ====================================================
%  VALIDATION PLOTS
% ====================================================
figure(2);

subplot(2,1,1)
plot(t_val, y1_val, 'b--', 'LineWidth', 1.3); hold on
plot(t_val, y1_val_pred, 'r', 'LineWidth', 2)
xlabel('Time'); ylabel('h_1');
title('FOPTD Validation – System 1 (h_1)');
legend('Validation data','Model prediction');
grid on

subplot(2,1,2)
plot(t_val, y2_val, 'b--', 'LineWidth', 1.3); hold on
plot(t_val, y2_val_pred, 'r', 'LineWidth', 2)
xlabel('Time'); ylabel('h_2');
title('FOPTD Validation – System 2 (h_2)');
legend('Validation data','Model prediction');
grid on
