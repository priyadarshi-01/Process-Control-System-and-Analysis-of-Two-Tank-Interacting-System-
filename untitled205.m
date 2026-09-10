clc; clear; close all;

%% ============================
%  Two-Tank Nonlinear System
%  Feedback Control Design
%  CV: h2 (level in tank 2)
%  MV: qi (inlet flow to tank 1)
%  DV: C2 (outlet valve coefficient of tank 2)
% =============================

%% -----------------------------
%  Physical Parameters
% ------------------------------
A1 = 3;    % ft^2
A2 = 5;    % ft^2
C1 = 8;    % cfm / ft^0.5
C2_nom = 12; % nominal C2 (disturbance variable)
qi_ss = 16;  % nominal inlet flow (cfm)

Ts = 1;          % sampling time (min)
Ns = 3000;       % total number of samples
k_vec = (1:Ns)'; % sample index (acts like "time" axis)

%% ============================
%  1. Steady-State of Nonlinear Plant
%     (analytical + numerical warm-up)
% ============================

% Analytical steady state (from mass balance)
h2_ss = (qi_ss / C2_nom)^2;
h1_ss = h2_ss + (qi_ss / C1)^2;

fprintf('Analytical steady-state:\n');
fprintf('  h1_ss = %.4f ft,  h2_ss = %.4f ft\n', h1_ss, h2_ss);

% (Optional) numerical check by integrating nonlinear model to steady state
h0 = [0; 0];                       % initial guess
tspan_ss = 0:Ts:500;              % long enough to reach SS
[~, h_ss_traj] = ode45(@(t,h) twoTank_nl(t,h,qi_ss,A1,A2,C1,C2_nom), ...
                       tspan_ss, h0);
h_ss_num = h_ss_traj(end,:);      % numerical SS
fprintf('Numerical steady-state (warm-up):\n');
fprintf('  h1_ss_num = %.4f ft,  h2_ss_num = %.4f ft\n', ...
        h_ss_num(1), h_ss_num(2));

% Use analytical SS as initial condition for closed-loop simulation
h1_0 = h1_ss;
h2_0 = h2_ss;

%% ============================
%  2. FOPTD Model & Direct Synthesis Tuning
%     (for h2 / qi)
% ============================

% FOPTD parameters for CV h2 (from your identification work)
Kp_proc = 0.1333;
Tau_P   = 0.9370;
Tau_D   = 0.0;

Tau_cl = Tau_P / 2;   % DO NOT change this if not allowed

% Direct Synthesis nominal gain
Kc_DS_nom = (1/Kp_proc) * (Tau_P / (Tau_cl + Tau_D));

% NEW: safety / robustness factor (0<alpha<1)
alpha = 0.3;                 % try 0.3, 0.4, 0.5 etc.

Kc_DS = alpha * Kc_DS_nom;   % reduced Kc without touching Tau_cl
Tau_I = Tau_P;               % or make integral slower: Tau_I = 1.5*Tau_P or 2*Tau_P


fprintf('\nDirect Synthesis controller settings (for h2/qi):\n');
fprintf('  Kc = %.4f\n', Kc_DS);
fprintf('  Tau_I = %.4f\n', Tau_I);

%% ============================
%  3. Closed-Loop Simulation
%     P controller and PI controller
% ============================

% ---------- Common Setpoint & Disturbance Profiles ----------
ysp_P  = zeros(Ns,1); % setpoint for P case
ysp_PI = zeros(Ns,1); % setpoint for PI case (same profile)

C2_prof = C2_nom * ones(Ns,1); % disturbance variable profile (C2)
% Servo: SP = h2_ss   for k=1..500
%        = h2_ss+5%   for k=501..3000
for k = 1:Ns
    if k <= 500
        ysp_P(k)  = h2_ss;
        ysp_PI(k) = h2_ss;
    else
        ysp_P(k)  = h2_ss * 1.05;  % +5% of SS
        ysp_PI(k) = h2_ss * 1.05;
    end
end

% Regulatory disturbance at k >= 2000: C2 increases by 5%
C2_prof(2000:end) = 1.05 * C2_nom;

% ============================
%  3A. Closed Loop with P Controller
% ============================

h1_P = zeros(Ns,1);  % tank 1 level trajectory
h2_P = zeros(Ns,1);  % tank 2 level trajectory (CV)
qi_P = zeros(Ns,1);  % manipulated variable trajectory
err_P = zeros(Ns,1); % error

% Initial conditions
h1_P(1) = h1_0;
h2_P(1) = h2_0;
qi_k    = qi_ss;

for k = 1:Ns
    % current setpoint and measurement
    y_sp = ysp_P(k);
    y_m  = h2_P(k);
    
    % error and P-control law
    err_P(k) = y_sp - y_m;
    qi_k = qi_ss + Kc_DS * err_P(k);   % P controller
    
    % (optional) saturation to avoid negative flow
    qi_k = max(qi_k, 0);
    
    % store MV
    qi_P(k) = qi_k;
    
    % current disturbance value (C2)
    C2_k = C2_prof(k);
    
    % plant simulation for one sampling interval [0 Ts]
    h_init = [h1_P(k); h2_P(k)];
    [~, h_traj] = ode45(@(t,h) twoTank_nl(t,h,qi_k,A1,A2,C1,C2_k), ...
                        [0 Ts], h_init);
    h_next = h_traj(end,:)';
    
    % store next state (unless last sample)
    if k < Ns
        h1_P(k+1) = h_next(1);
        h2_P(k+1) = h_next(2);
    end
end

%% ============================
%  3B. Closed Loop with PI Controller
% ============================

h1_PI = zeros(Ns,1);
h2_PI = zeros(Ns,1);
qi_PI = zeros(Ns,1);
err_PI = zeros(Ns,1);

% initial conditions
h1_PI(1) = h1_0;
h2_PI(1) = h2_0;
qi_k = qi_ss;
sum_err = 0;       % for integral action

for k = 1:Ns
    % current setpoint and measurement
    y_sp = ysp_PI(k);
    y_m  = h2_PI(k);
    
    % error
    err_PI(k) = y_sp - y_m;
    sum_err   = sum_err + err_PI(k);   % discrete integral (sum of errors)
    
    % PI control law (discrete form, same as CSTR example)
    qi_k = qi_ss + Kc_DS * ( err_PI(k) + (Ts/Tau_I)*sum_err );
    
    % saturation to avoid negative flow
    qi_k = max(qi_k, 0);
    
    % store MV
    qi_PI(k) = qi_k;
    
    % disturbance value
    C2_k = C2_prof(k);
    
    % plant simulation for one sampling interval
    h_init = [h1_PI(k); h2_PI(k)];
    [~, h_traj] = ode45(@(t,h) twoTank_nl(t,h,qi_k,A1,A2,C1,C2_k), ...
                        [0 Ts], h_init);
    h_next = h_traj(end,:)';
    
    if k < Ns
        h1_PI(k+1) = h_next(1);
        h2_PI(k+1) = h_next(2);
    end
end

%% ============================
%  4. Plots – Setpoint vs CV (h2)
% ============================

figure;

% ----- P Controller -----
subplot(2,1,1);
plot(k_vec, ysp_P, 'r', 'LineWidth', 1.8); hold on;
plot(k_vec, h2_P,  'b--',   'LineWidth', 1.8);
xlabel('Sample k');
ylabel('Level h_2 (ft)');
title('Closed-loop Response with P Controller');
legend('Setpoint','Controlled variable h_2','Location','Best');
grid on;

% ----- PI Controller -----
subplot(2,1,2);
plot(k_vec, ysp_PI, 'r', 'LineWidth', 1.8); hold on;
plot(k_vec, h2_PI,  'b--',   'LineWidth', 1.8);
xlabel('Sample k');
ylabel('Level h_2 (ft)');
title('Closed-loop Response with PI Controller');
legend('Setpoint','Controlled variable h_2','Location','Best');
grid on;

%% (Optional) Plot MVs and disturbance to see controller action
figure;
subplot(2,1,1);
plot(k_vec, qi_P, 'b--', 'LineWidth', 1.5);
xlabel('Sample k'); ylabel('q_i (P) [cfm]');
title('Manipulated variable q_i – P controller'); grid on;

subplot(2,1,2);
plot(k_vec, qi_PI, 'b--', 'LineWidth', 1.5);
xlabel('Sample k'); ylabel('q_i (PI) [cfm]');
title('Manipulated variable q_i – PI controller'); grid on;

%% ============================
%  END OF MAIN SCRIPT
% ============================


% ===========================================================
%  Nonlinear two-tank model (used by ode45)
%  A1, A2, C1, C2 are passed as parameters
% ===========================================================
function dhdt = twoTank_nl(~,h,qi,A1,A2,C1,C2)
    h1 = h(1); 
    h2 = h(2);

    % avoid negative arguments for sqrt
    delta = max(h1 - h2, 0);
    h2pos = max(h2, 0);

    q1 = C1 * sqrt(delta);   % flow from tank1 to tank2
    q2 = C2 * sqrt(h2pos);   % outlet of tank2

    dh1dt = (1/A1) * (qi - q1);
    dh2dt = (1/A2) * (q1 - q2);

    dhdt = [dh1dt; dh2dt];
end



