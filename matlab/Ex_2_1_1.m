clear
clc
close all


%% System invariant zeros
% Example from coursework
s = tf('s');
G = [1/s, 1/(s*(s+1));
    1/(s+1), 1/(s+1)];

[A, B, C, D] = ssdata(G);


[z, x_0, u_0] = system_zeros_square(A, B, C, D);

% Display results
disp('System Invariant Zeros (z):');
disp(z);
disp('Initial State Directions (x0):');
disp(x_0);
disp('Input Directions (u0):');
disp(u_0);


%% System transmission zeros
[Abar, Bbar, Cbar, Dbar] = minreal_custom(A, B, C, D);
[z, x_0, u_0] = system_zeros_square(Abar, Bbar, Cbar, Dbar);

% Display results
disp('System Transmission Zeros (z):');
disp(z);
disp('Initial State Directions (x0):');
disp(x_0); % Rotate the initial state vector to the initial space
disp('Input Directions (u0):');
disp(u_0); % The input vector direction doesn't change



%% Main Function
function [z, x_0, u_0] = system_zeros_square(A, B, C, D)

    % Construct Generaized Eigenvalue Problem
    [n, m] = size(B);
    M = [A, B; C, D];
    I_g = zeros(n+m);
    I_g(1:n,1:n) = eye(n);

    % Solve eigenvalue problem
    [V, D] = eig(M, I_g);
    eigenvalues = diag(D);
    
    % Remove infinite solutions
    finite_idx = isfinite(eigenvalues);
    z = eigenvalues(finite_idx);
    
    % Extract corresponding vectors
    full_vectors = V(:, finite_idx);
    x_0 = full_vectors(1:n, :);
    u_0 = full_vectors(n+1:end, :);
end


%% Helper functions
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


function Ct = ctrb(A, B)
    n = size(A, 1);
    
    % Initialize Ct with B (which is A^0 * B)
    Ct = B;
    
    % Current term tracks A^k * B to avoid re-calculating powers
    currentTerm = B;
    
    % Iterate from 1 to n-1 to append A*B, A^2*B, ...
    for k = 1:n-1
        currentTerm = A * currentTerm;
        Ct = [Ct, currentTerm];
    end
end


function Ob = obsv(A, C)
    n = size(A, 1);
    
    % Initialize Ob with C (which is C * A^0)
    Ob = C;
    
    % Current term tracks C * A^k
    currentTerm = C;
    
    % Iterate from 1 to n-1 to append C*A, C*A^2, ...
    for k = 1:n-1
        currentTerm = currentTerm * A;
        Ob = [Ob; currentTerm];
    end
end




