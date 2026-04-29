clear
clc
close all

%% 2.1.2 System zeros II
% System definition
s = tf('s');
G = 1/((s+1)*(s+2)*(s-1)) * ...
    [(s-1)*(s+2), 0, (s-1)^2;
    -(s+1)*(s+2), (s-1)*(s+1), (s-1)*(s+1)];

sys = ss(G); % This system is not minimal realization
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

[n, m] = size(B);

%% Define a search grid
% Select an outer box for the grid and a step size
real_min = -10; real_max = 10;
imag_min = -10; imag_max = 10;
step = 0.5; % Grid step

[Re, Im] = meshgrid(real_min:step:real_max, imag_min:step:imag_max);
S_grid = Re + 1i * Im;

fprintf('Checking %d points on the grid...\n', numel(S_grid));

%% Search inside the grid
threshold = 1e-5; % Below this value, the singular value is considered 0
candidates = [];

for k = 1:numel(S_grid)
    s = S_grid(k);
    
    % Rosenbrock Matrix
    P = [ (s*eye(n) - A), -B; C, D ];
    
    % Find the minimum singular value
    min_sv = min(svd(P));
    
    if min_sv < threshold
        candidates = [candidates; s];
    end
end

if isempty(candidates)
    error('No Singular Values in the grid below the threshold');
end

fprintf('%d zeros found.\n', length(candidates));


%% Directions and final print
disp('Final Results:');
fprintf('------------------------------------------------\n');


for i = 1:length(candidates)
    s_val = candidates(i);
    
    % Compute again the SVD to find the directions
    P = [(s_val*eye(n) - A), -B; C, D];
    [~, S, V] = svd(P); % We are only interested in the direction

    singular_vals = diag(S);
    r = sum(singular_vals > threshold);
    k = size(P, 2) - r;
    
    if k == 0
        disp('Warning: z0 does not appear to be a transmission zero.');
        X0_basis = []; U0_basis = [];
        return;
    end
    
    null_vec = V(:, end-k+1 : end);

    % Split the vector
    x0 = null_vec(1:n, :);
    u0 = null_vec(n+1:end, :);
    
    fprintf('Zero #%d: s = %.5f + %.5fi\n', i, real(s_val), imag(s_val));
    fprintf('x0:\n');
    disp(x0);
    fprintf('u0:\n');
    disp(u0);
    fprintf('------------------------------------------------\n');
end

% Final printing of the results
figure;
plot(real(S_grid(:)), imag(S_grid(:)), 'og','LineWidth', 1, 'MarkerSize', 1.5); hold on; % Grid
plot(real(candidates), imag(candidates), 'xr', 'LineWidth', 2, 'MarkerSize', 10); % Final zeros
legend('Grid', 'Zeros found');
grid on; title('Zeros search');
xlabel('Real'); ylabel('Imag');

