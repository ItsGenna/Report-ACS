clear
clc
close all

%% 2.2.3 Youla parametrization
s = tf('s');
G = (s-1)/(s^3+s^2-4*s-4);

Gss = minreal(ss(G));

A = Gss.a;
B = Gss.b;
C = Gss.c;
D = Gss.d;
% Calculate state - feedback and observer - gain matrices
F = -place(A, B, [-2 -5 -8]);
H = -transpose(place(A' , C' , [ -2 -5 -8]));
% Coprime factorisation
Gcf = tf(ss(A + B * F, [B, -H] , ...
    [F; C + D * F] , [1 , 0; D, 1]));
M = Gcf(1 ,1);
N = Gcf(2 ,1);
U = Gcf(1 ,2);
V = Gcf(2 ,2);

% Check if Bezout's identity holds
tf(minreal(V * M - U * N));
% Design controller
K = -simplify(zpk(minreal(U / V)));
step(feedback(G * K, 1)); grid on;

% Check internal stability
is_stable = ~any(real(find_poles(1/(1+G*K)))>0);

if is_stable
    fprintf("The controlled system is stable\n")
else
    fprintf("The controlled system is UNSTABLE\n")
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




