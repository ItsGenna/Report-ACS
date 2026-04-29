clear
clc
close all

%% State space model
A = [-2.2567e-02  -3.6617e+01  -1.8897e+01  -3.2090e+01   3.2509e+00  -7.6257e-01
       9.2572e-05  -1.8997e+00   9.8312e-01  -7.2562e-04  -1.7080e-01  -4.9652e-03
       1.2338e-02   1.1720e+01  -2.6316e+00   8.7582e-04  -3.1604e+01   2.2396e+01
       0            0   1.0000e+00            0            0            0 
       0            0            0            0  -3.0000e+01            0 
       0            0            0            0            0  -3.0000e+01];
B = [0     0 
      0     0 
      0     0 
      0     0 
     30     0 
      0    30];
C = [0     1     0     0     0     0 
      0     0     0     1     0     0];
D = [0     0 
      0     0];
sys = ss(A, B, C, D);

% Reduce model order
[hsv_stab, hsv_unstab] = hankelsv(sys, 'ncf', 'log'); % Small value correspond to states that do not contribute much to future state trajectory
sys_red = reduce(sys, 4, 'errortype', 'ncf'); % Reduce model to 4th-order


%% Get Smith form to find zeros
syms x
sys_tf = tf(sys_red);

% Find common denominator
common_den_poly = poly(sys_red.A); 
d_common = poly2sym(common_den_poly, x);


% Build numerator matrix
[num, den] = tfdata(sys_tf);
[rows, cols] = size(num);
N_mat = sym(zeros(rows, cols));
for i = 1:rows
    for j = 1:cols
        G_term = poly2sym(num{i,j}, x) / poly2sym(den{i,j}, x);
        N_mat(i,j) = simplify(G_term * d_common); 
    end
end

% Compute Smith Form
[U, V, Smith] = smithForm(N_mat, x);

% Find Smith-McMillan form
McMill = simplify(Smith / d_common);

% Get invariant zeros
diag_el = diag(McMill);
all_zeros = [];
all_poles = [];
for i=1:length(diag_el)
    [gamma, beta] = numden(diag_el(i));
    zeros_curr = solve(gamma == 0, x);
    all_zeros = [all_zeros; double(zeros_curr)];
    poles_curr = solve(beta==0, x);
    all_poles = [all_poles; double(poles_curr)];
end


% Display all zeros
disp("Zeros of the original system:");
disp(sort(tzero(sys)));

disp("Poles of the original system:");
disp(sort(pole(sys)));

disp("Zeros of the reduced system:");
disp(sort(tzero(sys_red)));

disp('Poles of the reduced system:');
disp(sort(pole(sys_red)));

disp('Zeros of the Smith-McMillan form:');
disp(sort(all_zeros));

disp('Poles of the Smith-McMillan form:');
disp(sort(all_poles));


% Plot Singular Values
figure;
sigma(sys, sys_red); 

hLines = findall(gcf, 'Type', 'line'); 
set(hLines, 'LineWidth', 1);
legend('Original System', 'Reduced System', 'Location', 'best');
title('Model Reduction Accuracy: Sigma Plot Comparison');
grid on;

%% Doubly-coprime factorization
A_red = sys_red.A;
B_red = sys_red.B;
C_red = sys_red.C;
D_red = sys_red.D;

F = - place(A_red, B_red, [-2, -5, -8, -10]);
H = - transpose(place(A_red' , C_red' , [-2, -5, -8, -10]));

A_F = A_red + B_red*F;
C_F = C_red + D_red*F;
A_H = A_red + H*C_red;
B_H = B_red + H*D_red;

M = ss(A_F, B_red, F, eye(size(B_red, 2)));
N = ss(A_F, B_red, C_F, D_red);
M_tilde = ss(A_H, H, C_red, eye(size(H, 2)));
N_tilde = ss(A_H, B_H, C_red, D_red);

X = ss(A_F, -H, C_F, eye(size(H, 2)));
Y = ss(A_F, -H, F, zeros(size(H, 2)));
X_tilde = ss(A_H, -B_H, F, eye(size(B_H, 2)));
Y_tilde = ss(A_H, -H, F, zeros(size(H, 2)));

% Check if Bezout's identity is satisfied
Bezout_Product = minreal([X_tilde, -Y_tilde; -N_tilde, M_tilde] * [M, Y; N, X]);
I_matrix = eye(size(Bezout_Product, 1));
error_norm = norm(Bezout_Product - I_matrix, inf);

disp(['Bezout Identity Error: ', num2str(error_norm)]);


%% Stabilizing controller
G = tf(sys_red);
Q = zeros(size(X));
K = - minreal((Y - M*Q) * inv(X - N*Q));

L = G*K;
S = inv(eye(size(L, 1)) + L);
T = minreal(L * S);

figure; sigma(S); title('Sensitivity (S)'); grid on;
figure; sigma(T); title('Complementary Sensitivity (T)'); grid on;
figure; step(feedback(G*K, eye(size(K, 1)))); grid on;


