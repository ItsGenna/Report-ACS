clear
clc
close all

%% 2.1.3 System perturbation
s = tf('s');
G = [7, 8; 6, 7] * [1/(s+1), 0; 0, 2/(s+2)] / [7, 8; 6, 7];
[n, m] = size(G);
I = eye(n);
K = I; % It should -I but with positive feedback, which is equivalent


%% a) Internal stability
tests = @(G,K) [inv(I+K*G), -K*inv(I+G*K), G*inv(I+K*G), inv(I+G*K)];

tests_int = tests(G, K);
isRH_inf = true;
for i = 1:length(tests_int)/n
    if ~is_RH_inf(tests_int(:, (i-1)*n+1:i*n))
        isRH_inf = false;
    end
end

if isRH_inf
    disp("The system is in RH_infinity");
else
    disp("The system is NOT in RH_infinity")
end

%% b) Stability margins
W1 = I;
W2 = I;

%% Additive uncertainty
M_add = W2 * (K / (I + G*K)) * W1;
[gamma_add, w_add] = hinfnorm(M_add);
margin_add = 1 / gamma_add;

% Display margin
fprintf('Additive Margin: %.4f at frequency %.4f rad/s\n', margin_add, w_add);

% Plot Singular Values
figure;
sigma(M_add);
grid on;
title('Singular Value Plot: Additive Uncertainty');

% Add a marker at the worst-case frequency
hold on;
plot(w_add, 20*log10(gamma_add), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
legend('Singular Values', ['Worst Case (\omega = ' num2str(w_add) ' rad/s)'], "Location","best");


%% Multiplicative uncertainty
M_mult = W2 * G * K / (I + G*K) * W1;
[gamma_mult, w_mult] = hinfnorm(M_mult);
margin_mult = 1 / gamma_mult;

% Diplay margin
fprintf('Multiplicative Margin: %.4f at frequency %.4f rad/s\n', margin_mult, w_mult);


% Plot Singular Values
figure;
sigma(M_mult);
grid on;
title('Singular Value Plot: Multiplicative Uncertainty');
hold on;

% Add a marker at the worst-case frequency
plot(w_mult, 20*log10(gamma_mult), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
legend('Singular Values', ['Worst Case (\omega = ' num2str(w_mult) ' rad/s)'], "Location","best");


%% c) Robust instability


%% Additive uncertainty
M_add_w0 = squeeze(freqresp(M_add, w_add));
[U_add, S_add, V_add] = svd(M_add_w0);

Delta_add = - gamma_add * V_add(:,1) * U_add(:,1)';
G_p_add = G + W1*Delta_add*W2;

% Get the matrices to test
tests_add = tests(G_p_add, K);

isRH_inf = true;
for i = 1:length(tests_add)/n
    if ~is_RH_inf(tests_add(:, (i-1)*n+1:i*n))
        isRH_inf = false;
        break
    end
end

if isRH_inf
    disp("The system with additive uncertainty is in RH_infinity");
else
    disp("The system with additive uncertainty is NOT in RH_infinity")
end


%% Multiplicative uncertainty
M_mult_w0 = squeeze(freqresp(M_mult, w_mult));
[U_mult, S_mult, V_mult] = svd(M_mult_w0);

Delta_mult = -gamma_mult * V_mult(:,1) * U_mult(:,1)';
G_p_mult = (I + W1*Delta_mult*W2) * G;
tests_mult = tests(G_p_mult, K);

isRH_inf = true;
for i = 1:length(tests_mult)/n
    if ~is_RH_inf(tests_mult(:, (i-1)*n+1:i*n))
        isRH_inf = false;
        break
    end
end
  
if isRH_inf
    disp("The system with multiplicative uncertainty is in RH_infinity");
else
    disp("The system with multiplicative uncertainty is NOT in RH_infinity")
end


%% Helper functions
function p = find_poles(H)
    syms s

    % First we need to convert from tf to symbolic function
    [num, den] = tfdata(H); 
    % Initialize an empty symbolic matrix
    [rows, cols] = size(H);
    H_sym = sym(zeros(rows, cols));
    
    % Loop through the matrix to convert each element
    for i = 1:rows
        for j = 1:cols
            % Convert coeff vector to symbolic polynomial and divide
            numerator_poly = poly2sym(num{i,j}, s);
            denominator_poly = poly2sym(den{i,j}, s);
            
            H_sym(i,j) = numerator_poly / denominator_poly;
        end
    end
    
    % Extract the numerator and denominator
    minor_list = [];
    for i=1:min(rows,cols)
        minor_list = [minor_list; compute_minors(H_sym, i)];
    end

    % Compute the pole polynomial
    [~, den_list] = numden(minor_list);
    pole_poly = prod(factor(den_list(1)));
    for k = 2:length(den_list)
        pole_poly = lcm(pole_poly, prod(factor(den_list(k))));
    end

    p = double(solve(pole_poly, s));
end


function det_list = compute_minors(G, ord)
    % Computes all the minors of G of order ord
    
    % Get system dimensions
    [p,m] = size(G);

    % Handle special cases
    if p<ord || m<ord
        det_list = [];
        return
    end
    
    if (p==ord && m==ord)
        det_list = simplify(prod(factor(det(G))));
        return
    end

    
    % Get all combinations of column and row indices to KEEP
    col_combos = nchoosek(1:m, ord);
    row_combos = nchoosek(1:p, ord);
    

    % Loop through combinations to compute determinants
    num_row_sets = size(row_combos, 1);
    num_col_sets = size(col_combos, 1);
    det_list = sym(zeros(num_row_sets * num_col_sets, 1));

    counter = 1;
    for i = 1:num_row_sets
        for j = 1:num_col_sets
            % Extract the specific indices for this iteration
            r_idx = row_combos(i, :);
            c_idx = col_combos(j, :);
            
            % Extract the submatrix
            sub_G = G(r_idx, c_idx);
            
            % Compute determinant and store it
            det_list(counter) = simplify(prod(factor(det(sub_G))));
            counter = counter + 1;
        end
    end
end


function answer = is_RH_inf(H)
    answer = true;
    is_Stable = true;
    if any(real(find_poles(H))>=0)
        is_Stable = false;
    end

    % Check for Properness
    isProper = isproper(H);

    % Check for Rationality
    isRational = ~hasdelay(H);

    if ~(is_Stable && isProper && isRational)
        answer = false;
    end
end




