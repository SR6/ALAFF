function T = Spectral_Decomposition_Lambda( T )
% Returns Lambda such that T = Q Lambda Q' is the
%  Spectral Decomposition of T

    tol = 1e-15;  % Deflation threshold: used 1e-10 to 1e-15 with not much difference
    m = size(T, 1);

    % Work with a copy of T to avoid overwriting upper triangle
    while m > 2
        % Check for deflation from bottom up (smallest eigenvalue to
        % isolate first)
        while m > 2 && abs(T(m, m-1)) < tol
            m = m - 1;  % Deflate if appropriate
        end

        % Apply Francis Step to the leading m x m submatrix
        T(1:m, 1:m) = Francis_Step(T(1:m, 1:m));
    end

    % Handle final 2x2 block using eig
    A = T(1:2, 1:2);  % Extract the final 2x2 block
    lambda = eig(A);  % Compute eigenvalues
    T(1,1) = lambda(1);
    T(2,2) = lambda(2);
    T(2,1) = 0;  % Zero out subdiagonal

    % Final cleanup: zero out any tiny subdiagonals left
    for i = 2:size(T,1)
        if abs(T(i,i-1)) < tol
            T(i,i-1) = 0;
        end
    end
end
