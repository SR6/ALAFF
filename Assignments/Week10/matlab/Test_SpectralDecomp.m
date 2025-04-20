format long

rng(0);

m = 6;

% Create a random symmetric tridiagonal matrix
T = rand(m, m);
T = triu(T, -1);
T = tril(T) + tril(T, -1)';

T_original = T;

% Call spectral decomposition
T_diag = Spectral_Decomposition_Lambda(T);

% Display the final matrix to inspect convergence
disp('Final matrix T_diag from Spectral_Decomposition_Lambda:');
disp(T_diag);

% Print subdiagonal entries to check convergence toward zero
fprintf('\nSubdiagonal entries T_diag(i+1,i):\n');
for i = 1:m-1
    fprintf('T_diag(%d,%d) = %.3e\n', i+1, i, T_diag(i+1, i));
end

% Extract eigenvalues from result (diagonal)
lambda_custom = sort(diag(T_diag));
lambda_true = sort(eig(T_original));

% Display results
disp('Eigenvalues from Spectral_Decomposition_Lambda:');
disp(lambda_custom');

disp('Eigenvalues from eig():');
disp(lambda_true');

% Error norm
fprintf('L2 norm of eigenvalue error: %.3e\n', norm(lambda_custom - lambda_true));

% Optional: plot
figure;
stem(lambda_custom, 'filled'); hold on;
stem(lambda_true, 'r');
legend('Spectral Decomposition', 'eig()');
title('Eigenvalue Comparison');
xlabel('Index'); ylabel('Eigenvalue'); grid on;
