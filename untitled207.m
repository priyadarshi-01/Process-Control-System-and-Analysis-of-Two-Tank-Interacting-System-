clc; clear; close all;

%% ============================================
% CLOSED-LOOP STABILITY ANALYSIS – P CONTROLLER
% Process: h2(s) / qi(s)
% ============================================

% Identified process parameters
Kp  = 0.1333;
tau = 0.937;

%% ============================================
% (a) ROUTH–HURWITZ STABILITY ANALYSIS
% ============================================

syms Kc s real

% Closed-loop characteristic equation:
% 1 + Kc*Kp/(tau*s + 1) = 0
char_eq = tau*s + (1 + Kc*Kp);

disp('Closed-loop characteristic equation:')
pretty(char_eq)

% Stability condition from Routh-Hurwitz
Kc_min = -1 / Kp;

fprintf('\n======= ROUTH–HURWITZ RESULT =======\n');
fprintf('Closed-loop system is stable for:\n');
fprintf('Kc > %.4f\n', Kc_min);
fprintf('Hence, system is stable for all positive Kc.\n');

%% ============================================
% (b) ROOT LOCUS ANALYSIS
% ============================================

% Open-loop transfer function (without gain)
Gp = tf(Kp, [tau 1]);

figure;
rlocus(Gp)     
grid on;
title('Root Locus of Closed-Loop System with P Controller')
xlabel('Real Axis')
ylabel('Imaginary Axis')

%% ============================================
% Interpretation (Command Window Output)
% ============================================

fprintf('\n======= ROOT LOCUS INTERPRETATION =======\n');
fprintf('• Open-loop has one pole at s = %.4f\n', -1/tau);
fprintf('• No zeros present\n');
fprintf('• Root locus lies entirely on negative real axis\n');
fprintf('• Increasing Kc moves pole further left\n');
fprintf('• Closed-loop system is stable for all Kc > 0\n');
