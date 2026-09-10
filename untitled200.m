clc; clear; close all;
%% Assignment 01 (modified for up-then-down steps +10%→-10% and +20%→-20%)
%% Step 0: Define Parameters ---
A1 = 3; A2 = 5;
C1 = 8; C2 = 12;
qi_ss = 16; % steady-state inlet flow (cfm)
Ts = 1;               % sampling time (min)
tspan = 0:Ts:200;     % total simulation time (200 min)

%% Step 1: Steady-State Calculation ---
h2_ss = (qi_ss/C2)^2;
h1_ss = h2_ss + (qi_ss/C1)^2;
fprintf('Steady state levels:\n  h1_ss = %.3f ft, h2_ss = %.3f ft\n', h1_ss, h2_ss);

%% Step 2: Define Nonlinear Model ---
twoTank = @(t,h,qi) [
    (1/A1)*(qi - C1*sqrt(max(h(1)-h(2),0)));
    (1/A2)*(C1*sqrt(max(h(1)-h(2),0)) - C2*sqrt(max(h(2),0)))
];

%% Step 3: Transient to Steady-State (optional check) ---
h0 = [0; 0];
tspan_trans = 0:Ts:200;
[~, h_trans] = ode45(@(t,h) twoTank(t,h,qi_ss), tspan_trans, h0);

figure;
plot(tspan_trans, h_trans(:,1), 'b', 'LineWidth', 1.5); hold on;
plot(tspan_trans, h_trans(:,2), 'r', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Liquid Level (ft)');
title('Transient to Steady-State Profile for Two-Tank System');
legend('h_1(t)','h_2(t)'); grid on;

%% Step 4: Nonlinear responses for up-then-down steps (with noise)
step_percent = [0.10 0.20]; 
colors = {'r','b'};
h_nonlinear_all = cell(length(step_percent),1);

% Define step timing:
t_pre  = 0:Ts:100;            % steady (pre)
t_up   = 100:Ts:150;          % up-step hold (50 min)
t_down = 150:Ts:200;          % down-step hold (50 min)

figure;
subplot(2,1,1); hold on;
subplot(2,1,2); hold on;

for k = 1:length(step_percent)
    X = step_percent(k);
    %% Pre-step (0-100)
    [~, h_pre] = ode45(@(t,h) twoTank(t,h,qi_ss), t_pre, [0;0]);
    
    %% Up-step (100-150): +X%
    qi_up = qi_ss * (1 + X);
    [~, h_up] = ode45(@(t,h) twoTank(t,h,qi_up), t_up, h_pre(end,:));
    
    %% Down-step (150-200): -X% (relative to steady-state)
    qi_down = qi_ss * (1 - X);
    [~, h_down] = ode45(@(t,h) twoTank(t,h,qi_down), t_down, h_up(end,:));
    
    %% Combine
    t_full = [t_pre, t_up(2:end), t_down(2:end)];
    h_full = [h_pre; h_up(2:end,:); h_down(2:end,:)];
    
    %% Add noise (small)
    h1_noisy = h_full(:,1) + 0.01*h1_ss*randn(size(h_full(:,1)));
    h2_noisy = h_full(:,2) + 0.01*h2_ss*randn(size(h_full(:,2)));
    
    h_nonlinear_all{k} = [t_full(:), h1_noisy(:), h2_noisy(:)];
    
    %% Plot noisy nonlinear
    subplot(2,1,1);
    plot(t_full, h1_noisy, colors{k}, 'LineWidth', 1.2);
    ylabel('h_1 (ft)');
    title('Nonlinear Two-Tank Response (with noise) - h_1');
    grid on;
    
    subplot(2,1,2);
    plot(t_full, h2_noisy, colors{k}, 'LineWidth', 1.2);
    ylabel('h_2 (ft)'); xlabel('Time (min)');
    title('Nonlinear Two-Tank Response (with noise) - h_2');
    grid on;
end

subplot(2,1,1); legend('+10%→-10%','+20%→-20%');
subplot(2,1,2); legend('+10%→-10%','+20%→-20%');

%% Step 5: Symbolic Linearization using Jacobians including q2 ---
syms h1 h2 qi q2 positive

q1 = C1*sqrt(h1 - h2);
f1 = (1/A1)*(qi - q1);
f2 = (1/A2)*(q1 - q2);
f = [f1; f2];

J_f_states_inputs = jacobian(f, [h1 h2 qi q2]);
q2_expr = C2*sqrt(h2);
J_f_states_inputs = subs(J_f_states_inputs, q2, q2_expr);

q_vec = [q1; q2_expr];
J_q_states = jacobian(q_vec, [h1 h2]);

q2_ss = C2*sqrt(h2_ss);
subs_vars = [h1 h2 qi];
subs_vals = [h1_ss h2_ss qi_ss];

A_extended = double(subs(J_f_states_inputs, subs_vars, subs_vals));
Jq_states  = double(subs(J_q_states, subs_vars(1:2), subs_vals(1:2)));

fprintf('\nJacobian A_extended = ∂f/∂[h1,h2,qi,q2]:\n'); disp(A_extended);
fprintf('Jacobian of flows w.r.t states:\n'); disp(Jq_states);

%% Step 6: State-space model creation ---
C1_lin = [1 0];
C2_lin = [0 1];
D = 0;

A = A_extended(:,1:2);
B = A_extended(:,3);

sys_lin_h1 = ss(A, B, C1_lin, D);
sys_lin_h2 = ss(A, B, C2_lin, D);

%% Step 7: Linear vs Nonlinear Comparison (up→down) ---
figure;
for var = 1:2
    subplot(2,1,var); hold on;
    for k = 1:length(step_percent)
        X = step_percent(k);
        % deviation input u(t)
        u = zeros(length(tspan),1);
        u(tspan>=100 & tspan<150) = X*qi_ss;
        u(tspan>=150)             = -X*qi_ss;
        
        if var == 1
            [y_lin, t_lin] = lsim(sys_lin_h1, u, tspan);
            h_lin = y_lin;        
            h_non = h_nonlinear_all{k}(:,2) - h1_ss;
            ylabel_text = 'h_1 deviation (ft)';
            title_text  = 'Linear vs Nonlinear Response (Deviation) – h_1';
        else
            [y_lin, t_lin] = lsim(sys_lin_h2, u, tspan);
            h_lin = y_lin;
            h_non = h_nonlinear_all{k}(:,3) - h2_ss;
            ylabel_text = 'h_2 deviation (ft)';
            title_text  = 'Linear vs Nonlinear Response (Deviation) – h_2';
        end
        
        plot(t_lin, h_lin, colors{k}, 'LineWidth', 1.5);
        plot(h_nonlinear_all{k}(:,1), h_non, [colors{k} '--'], 'LineWidth', 1.2);
    end
    xlabel('Time (min)'); ylabel(ylabel_text);
    title(title_text);
    legend('Linear +10%(up)','Linear +20%(up)', ...
           'Nonlinear +10%(up)','Nonlinear +20%(up)');
    grid on;
end

%% Step 8: Transfer Function & Stability ---
[num, den] = ss2tf(A,B,C2_lin,D);   % h2/qi
G_tf = tf(num,den);
disp('---------------------------------------------');
disp('Transfer Function G(s)/Q_i(s):');
G_tf

poles = eig(A);
disp('Poles of A:'); disp(poles);
if all(real(poles)<0)
    disp('System is stable.');
else
    disp('System is unstable.');
end

%% Step 8B — Transfer Function G(s)/Q_2(s) ---
C_q2 = [0  C2/(2*sqrt(h2_ss))];
D_q2 = 0;

[num_q2, den_q2] = ss2tf(A, B, C_q2, D_q2);
G_q2 = tf(num_q2, den_q2);

disp('---------------------------------------------');
disp('Transfer Function G(s)/Q_2(s): (Outlet flow / Inlet flow)');
G_q2

poles_q2 = pole(G_q2);
disp('Poles of G(s)/Q_2(s):'); disp(poles_q2);
if all(real(poles_q2)<0)
    disp('Outlet-flow transfer function is stable.');
else
    disp('Outlet-flow transfer function is unstable.');
end

%% Step 9: Compute Kp and Tau (ignoring integrator pole at s=0)
systems = {G_tf, 'G_Qi'; ...
           G_q2, 'G_Q2'};

for i = 1:size(systems,1)
    sys  = systems{i,1};
    name = systems{i,2};

    % Convert to zero-pole-gain
    [z, p, k] = zpkdata(sys,'v');

    % Remove the integrator pole at s = 0
    tol = 1e-6;
    p_no_int = p(abs(p) > tol);   % remaining (stable) pole

    % First-order part (no integrator)
    sys_fo = zpk(z, p_no_int, k);

    % Process gain Kp and time constant tau
    Kp  = dcgain(sys_fo);         % steady-state gain of first-order part
    tau = -1/real(p_no_int(1));   % tau = -1/p  (p is negative)

    fprintf('%s: Kp = %.4f, Tau = %.4f\n', name, Kp, tau);
end
