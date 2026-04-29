clear
clc

%% Find the uncertainty representation
% [a, b, c, d] = uncertainty();

% fprintf("a: %.3f\nb: %.3f\nc: %.3f\nd: %.3f\n", a, b, c, d);


%% Find Uncertainty matrices with Upper LFT
% a = 2.5; b = -0.7; c = 1; d = -0.6;
% 
% M = upper_lft(a, b, c, d);
% disp(M)


%% Find a minimal realization of the state space representation of W1 and W2
s = tf('s');
W1 = 2 / (s+0.1127);
W2 = 0.3 * (0.7*s+1) / (0.07*s+1);

sys1 = minreal(ss(W1));
sys2 = minreal(ss(W2));

[Aw1, Bw1, Cw1, Dw1] = ssdata(sys1);
[Aw2, Bw2, Cw2, Dw2] = ssdata(sys2);


%% Uncertainty channels with Control Systems Toolbox
m = ureal('m', 10, 'Range', [5, 18]);
c = ureal('c', 10, 'Range', [6, 14]);
k = ureal('k', 1, 'Range', [0.8, 1.8]);
h = ureal('h', 2.5, 'Range', [2, 4.5]);

A = [0, 1; -k/m, -c/m];
B = [0; 1/m];
C = [h, 0];
D = 0;

G_p = uss(A, B, C, D);

%% LFT Automatic
A_N = [0, 1, 0, 0;
    -k/m, -c/m, 0, 0;
    Bw1*h, 0, Aw1, 0;
    0, 0, 0, Aw2];

B_N = [0, 0;
    0, 1/m;
    Bw1, 0;
    0, Bw2];

C_N = [Dw1*h, 0, Cw1, 0;
    0, 0, 0, Cw2;
    -h, 0, 0, 0];

D_N = [Dw1, 0;
    0, Dw2;
    -1, 0];

Plant_aut = ss(A_N, B_N, C_N, D_N);
[P_aut, Delta_aut] = lftdata(Plant_aut);


%% LFT Manual
M_m = upper_lft(1, -0.231, 10, 3.846);
M_k = upper_lft(1, -0.28, 1, -0.6);
M_h = upper_lft(2.5, -0.7, 1, -0.6);

A_P = [0, 1, 0, 0;
    -M_m(2, 2)*M_k(2, 2), -M_m(2, 2)*c.NominalValue, 0, 0;
    Bw1*M_h(2, 2), 0, Aw1, 0;
    0, 0, 0, Aw2];

B_P = [0, 0, 0, 0, 0, 0;
    M_m(2, 1), -M_m(2, 2), -M_k(2, 1)*M_m(2, 2), 0, 0, M_m(2, 2);
    0, 0, 0, Bw1*M_h(2, 1), Bw1, 0;
    0, 0, 0, 0, 0, Bw2];

C_P = [-M_m(1, 2)*M_k(2, 2), -M_m(1, 2)*c.NominalValue, 0, 0;
    0, 0.4*c.NominalValue, 0, 0;
    M_k(1, 2), 0, 0, 0;
    M_h(1, 2), 0, 0, 0;
    Dw1*M_h(2, 2), 0, Cw1, 0;
    0, 0, 0, Cw2;
    -M_h(2, 2), 0, 0, 0];

D_P = [M_m(1, 1), -M_m(1, 2), -M_m(1, 2)*M_k(2, 1), 0, 0, M_m(1, 2);
    0, 0, 0, 0, 0, 0;
    0, 0, M_k(1, 1), 0, 0, 0;
    0, 0, 0, M_h(1, 1), 0, 0;
    0, 0, 0, Dw1*M_h(2, 1), Dw1, 0;
    0, 0, 0, 0, 0, Dw2;
    0, 0, 0, -M_h(2, 1), -1, 0];

P = ss(A_P, B_P, C_P, D_P);



%% Design a controller (and Nominal Performance)
G_nom = tf(minreal(G_p.NominalValue));
K = 0.3 * (s/0.88729 + 1) * (s/0.1127 + 1) * (s/48 + 1)/ (s*(s/2 + 1)^2);

% Loop Transfer Function
L = G_nom * K;

% Sensitivity Function
S = 1 / (1 + L); % (Transfer function from Reference to Error)

% Control Sensitivity
T_ru = K * S; % (Transfer function from Reference r to Control Signal u)


% Plot the Bode Diagrams
figure('Units', 'normalized', 'Position', [0.3, 0.1, 0.45, 0.8]);

subplot(2,1,1);
bodemag(S, 'b', 1/W1, 'r--');
grid on;
title('Weighted Sensitivity (W1 * S)');
legend('S', 'Target (1/W1)', 'Location', 'Best');

subplot(2,1,2);
bodemag(T_ru, 'b', 1/W2, 'r--');
grid on;
title('Weighted Control Sensitivity (W2 * K * S)');
legend('K*S', 'Target (1/W2)', 'Location', 'Best');


%% M-Delta structure
M = lft(P, K);
fprintf("The norm of M22 computed manually is: %f\n", hinfnorm(minreal(M(5:6, 5))));

M_aut = lft(P_aut, K);
fprintf("The norm of M22 computed automatically is: %f\n", hinfnorm(M_aut(7:8, 7)));

%% Compare the two models
Delta_m = ureal('Delta_m', 0, 'Range', [-1, 1]);
Delta_c = ureal('Delta_c', 0, 'Range', [-1, 1]);
Delta_k = ureal('Delta_k', 0, 'Range', [-1, 1]);
Delta_h = ureal('Delta_h', 0, 'Range', [-1, 1]);
Delta = [Delta_m, 0, 0, 0;
        0, Delta_c, 0, 0;
        0, 0, Delta_k, 0;
        0, 0, 0, Delta_h];
F = lft(Delta, M);


F_aut = lft(Delta_aut, M_aut);

figure();
sigma(F.NominalValue);
hold on; grid on;
sigma(F_aut.NominalValue);

% Compare transfer functions
tf_manual = tf(F);
tf_aut = tf(F_aut);

figure();
bode(tf_manual(1));
hold on; grid on;
bode(tf_aut(1));

figure();
bode(tf_manual(2));
hold on; grid on;
bode(tf_aut(2));


%% Robust stability
fprintf("The norm of M11 computed manually is: %f\n", hinfnorm(minreal(M(1:4, 1:4))));
fprintf("The norm of M11 computed automatically is: %f\n", hinfnorm(minreal(M_aut(1:6, 1:6))));


%% Less conservative procedure
% Manual system
Nfr = frd(minreal(M(1:4, 1:4)), logspace(-3, 3, 100));
[bounds, muInfo] = mussv(Nfr, [-1, 1; -1, 1; -1, 1; -1, 1]);
mu = bounds(1);

if all(mu.ResponseData<1)
    disp("The system is robustly stable with structured singular value");
else
    disp("The system may NOT be robustly stable with structured singulaar value");
end


% Control System Toolbox system
Nfr_aut = frd(minreal(M_aut(1:4, 1:4)), logspace(-3, 3, 100));
[bounds_aut, muInfo_aut] = mussv(Nfr_aut, [-1, 0; -1, 0; -1, 0; -1, 0]);
mu_aut = bounds_aut(1);

if all(mu_aut.ResponseData<1)
    disp("The system is robustly stable with structured singular value");
else
    disp("The system may NOT be robustly stable with structured singular value");
end



%% Pade approximation
tau = ureal('tau', 0.1, 'perc', [-30, 30]);
pad = (1-tau/2*s) / (1+tau/2*s);
G = h/(m*s^2 + c*s + k) * pad;
S = uss(1/(1 + G*K));
G_delay = uss([W1*S; ...
        W2*K*S]);

[M_delay, Delta_delay] = lftdata(G_delay);
n_delay = length(Delta_delay.NominalValue);
M_11 = minreal(M_delay(1:n_delay, 1:n_delay));

if hinfnorm(M_11)<1
    disp("The system with time delay is robustly stable with small gain theorem");
else
    disp("The system with time delay may NOT be robustly stable with small gain theorem");
end


% Less conservative procedure
Nfr = frd(minreal(M_11), logspace(-3, 3, 100));
[bounds_delay, muInfo_delay]= mussv(Nfr, [-1, 1; -1, 1; -1, 1; -1, 1;-1, 1;-1, 1;-1, 1;-1, 1;-1, 1;-1, 1;-1, 1;-1, 1;]);
mu_delay = bounds_delay(1);

if all(mu_delay.ResponseData<1)
    disp("The system with time delay is robustly stable with structured singular value");
else
    disp("The system with time delay may NOT be robustly stable with structured singular value");
end


%% Helper functions
function [a, b, c, d] = uncertainty(a_max, a_min, a_nom)
    syms a b c d
    [a, b, c, d] = solve((a+b)/(c+d)==a_max, (a-b)/(c-d)==a_min, (a/c)==a_nom, c==1);

    a = eval(a); b = eval(b); c = eval(c); d = eval(d);
end

function M = upper_lft(a, b, c, d)
    M = [-d/c, 1/c; (-a*d/c+b), a/c];
end
