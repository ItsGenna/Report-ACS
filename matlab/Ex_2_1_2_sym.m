clear
clc
close all

%% 2.1.2 System zeros II
s = tf('s');
G = 1/((s+1)*(s+2)*(s-1)) * ...
    [(s-1)*(s+2), 0, (s-1)^2;
    -(s+1)*(s+2), (s-1)*(s+1), (s-1)*(s+1)];

sys = ss(G);
[z, x_0, u_0] = zeros_non_square(sys.A, sys.B, sys.C, sys.D);


% Display results
disp('System Zeros (z):');
disp(z);
disp('Initial State Directions (x0):');
disp(x_0);
disp('Input Directions (u0):');
disp(u_0);

[Abar, Bbar, Cbar, Dbar] = minreal_custom(sys.A, sys.B, sys.C, sys.D);




%% Main function
function [z, x_0, u_0] = zeros_non_square(A, B, C, D)
    syms s

    % Get system dimensions
    [n, m] = size(B);
    p = size(C, 1);
    I = eye(n);

    % Compute symbolic transfer function matrix
    G = C*((s*I - A)\ B) + D;
    
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
        % Determine Geometric Multiplicity (basis of null space)
        % Z_block = null(P); 
        [X0_basis, U0_basis] = get_multizero_directions(A, B, C, D, z_unique(i));

        % Store results
        x_0 = X0_basis; %[x_0, Z_block(1:n, :)];
        u_0 = U0_basis; %[u_0, Z_block(n+1:end, :)];
    end
end

%% Helper functions
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
    r = sum(singular_vals > tol);
    k = size(P, 2) - r;
    
    if k == 0
        disp('Warning: z0 does not appear to be a transmission zero.');
        X0_basis = []; U0_basis = [];
        return;
    end
    
    fprintf('Zero at s = %.2f has total null dimension k = %d\n', real(z0), k);

    % Extract the Subspace (Last k columns of V)
    V_null = V(:, end-k+1 : end);
    
    % Separate state (x) and input (u)
    X0_basis = V_null(1:n, :);
    U0_basis = V_null(n+1:end, :);
end


function [Abar, Bbar, Cbar, Dbar] = minreal_custom(A, B, C, D)
    % Returns a minimal realization of a state-space system
    % Get system sizes
    [n, ~] = size(B);

    % Compute controllability and observability matrices
    Ct = ctrb(A, B);
    Ob = obsv(A, C);

    % Calculate SVD of the controllability matrix
    [U, S, ~] = svd(Ct);

    % Get the rank of the controllability matrix
    rc = length(find(S>1e-12));
    Uc = U(1:n, 1:rc);
    Sc = S(1:rc, 1:rc);

    % Compute SVD again
    clear U S V
    [~, S, V] = svd(Ob * Uc * sqrtm(Sc));

    % Calculate rank
    rco = length(find(S>1e-12));
    Sco = S(1:rco, 1:rco);
    Vco = V(1:rc, 1:rco);

    % Setup transformation matrices
    R1 = (Uc / sqrtm(Sc)) * Vco * sqrtm(Sco);
    T1 = Uc * sqrtm(Sc) * (Vco / Sco);

    % Transform the system
    Abar = R1' * A * T1;
    Bbar = R1' * B;
    Cbar = C * T1;
    Dbar = D;

end

