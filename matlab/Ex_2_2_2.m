clear
clc
close all

%% 
A = [2, 0; 0, -2];
B = [1, 0; 0, 1];
C = [1, -1];
D = [1, 1];


[z, x_0, u_0] = zeros_non_square(A, B, C, D);

disp("Zeros found:");
disp(z);

%% Helper functions
function [z, x_0, u_0] = zeros_non_square(A, B, C, D)
    syms s

    % Get system dimensions
    [n, m] = size(B);
    p = size(C, 1);
    I = eye(n);

    % Compute symbolic transfer function matrix
    G = simplifyFraction(C*((s*I - A)\ B) + D);
    
    % Extract the numerator and denominator
    minor_list = [];
    for i=1:min(m,p)
        minor_list = [minor_list; compute_minors(G, i)];
    end

    % Compute the pole polynomial
    [~, den_list] = numden(minor_list);
    pole_poly = prod(factor(den_list(1)));
    for k = 2:length(den_list)
        pole_poly = lcm(pole_poly, prod(factor(den_list(k))));
    end
    
    % Compute the numerators as if they had the pole polynomial at the
    % denominator
    zero_minor_list = compute_minors(G, rank(G));
    num_list = simplify(zero_minor_list * pole_poly);
    
    % Find the zero polynomial
    zero_poly = num_list(1); 
    for k = 2:length(num_list)
        zero_poly = gcd(zero_poly, num_list(k));
    end
    
    % Find the roots of the zero polynomial to get the zeros
    z = double(solve(zero_poly, s));
    z_unique = unique(z);

    x_0 = [];
    u_0 = [];


    for i = 1:length(z_unique)
        current_z = z_unique(i);
        
        % Define Rosenbrock Matrix
        P = [current_z*I - A, -B;
             C,             D];
        
        % Determine Geometric Multiplicity (basis of null space)
        % Z_block = null(P); 
        [X0_basis, U0_basis] = get_multizero_directions(A, B, C, D, z_unique(i));

        % Store results
        x_0 = X0_basis; %[x_0, Z_block(1:n, :)];
        u_0 = U0_basis; %[u_0, Z_block(n+1:end, :)];
    end


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


function [X0_basis, U0_basis] = get_multizero_directions(A, B, C, D, z0)
    % Inputs: A, B, C, D matrices and the specific zero location z0
    
    n = size(A, 1);
    
    % Form Rosenbrock Matrix
    P = [ (z0*eye(n) - A),  B ; 
          -C,               D ];
      
    % Perform SVD
    [~, S, V] = svd(P);
    
    % Determine Rank Deficiency
    % Find how many singular values are effectively zero
    singular_vals = diag(S);
    tol = 1e-5; % Tolerance for numerical zero
    
    % The indices of the 'zero' singular values (at the end of the list)
    zero_indices = find(singular_vals < tol);
    k = length(zero_indices);
    
    if k == 0
        disp('Warning: z0 does not appear to be a transmission zero.');
        X0_basis = []; U0_basis = [];
        return;
    end
    
    fprintf('Zero at s = %.2f has geometric multiplicity k = %d\n', real(z0), k);

    % Extract the Subspace (Last k columns of V)
    V_null = V(:, end-k+1 : end);
    
    % Separate state (x) and input (u)
    X0_basis = V_null(1:n, :);
    U0_basis = V_null(n+1:end, :);
    
    % Orthonormalize the bases
    X0_basis = orth(X0_basis);
    U0_basis = orth(U0_basis);
end

