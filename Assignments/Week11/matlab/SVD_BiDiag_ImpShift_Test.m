function [S, U, V] = SVD_BiDiag_ImpShift_Test(B)
    % Computes the SVD of a bidiagonal matrix B using implicit shift algorithm
    % Input: B - upper bidiagonal matrix
    % Output: S - vector of singular values (sorted in descending order)
    %         U, V - orthogonal matrices such that B = U * diag(S) * V'

    disp("in right file")
    [n, m] = size(B);
    min_dim = min(n, m);
    U = eye(n);
    V = eye(m);
    Bk = B;
    max_iter = 50000;% * max(n, m);
    tol = 1e-15;

    % disp("original B");
    % B
    for iter = 1:max_iter
        % Check for convergence
        for i = 1:min_dim - 1
            if abs(Bk(i, i + 1)) <= tol * (abs(Bk(i, i)) + abs(Bk(i + 1, i + 1)))
                Bk(i, i + 1) = 0;
            end
        end
        off_diag = abs(diag(Bk, 1));
        if isempty(off_diag) || max(off_diag) < tol * norm(Bk, 'fro') %remove max(off_diag)<tol?
            disp("reached convergence")
            break;
        end

        % Find the unreduced part of the bidiagonal matrix
        q = min_dim;
        while q > 1 && abs(Bk(q - 1, q)) <= tol
            q = q - 1;
            q
        end
        if q <= 1
            continue; % Nothing to do
        end
        p = 1;
        while p < q - 1 && abs(Bk(p, p + 1)) <= tol
            p = p + 1;
            p
        end

        % Wilkinson shift (for B^T * B)
        mu = 0;
        if q > 1
            s_qq = Bk(q, q)^2;
            s_qq_minus_1 = Bk(q - 1, q - 1)^2;
            f_qq_minus_1 = Bk(q - 1, q)^2;
            delta = (s_qq_minus_1 - s_qq) / (2 * f_qq_minus_1);
            mu = s_qq + delta - sign(delta) * sqrt(delta^2 + 1) * f_qq_minus_1;
        else
            mu = Bk(1, 1)^2;
        end

        % Initial Givens rotation to create the bulge
        x = Bk(p, p)^2 - mu;
        y = Bk(p, p) * Bk(p, p + 1);
        G_right_0 = Givens_rotation([x; y]);
        Bk(p:min(p + 1, m), p:p + 1) = Bk(p:min(p + 1, m), p:p + 1) * G_right_0;
        V(:, p:p+1) = V(:, p:p+1) * G_right_0;

        % disp("after initial given's")
        % Bk
        % Chase the bulge
        for k = p:q - 1
            % Left rotation to eliminate the subdiagonal (if it exists)
            if k < n %- 1
                G_left = Givens_rotation([Bk(k, k); Bk(k + 1, k)]);
                Bk(k:k + 1, k:min(m, k + 2)) = G_left' * Bk(k:k + 1, k:min(m, k + 2));
                U(:, k:k + 1) = U(:, k:k + 1) * G_left;
            end

            % disp("after left elim sub diagonal")
            % Bk
            % Right rotation to eliminate the created superdiagonal
            if k < q && k + 2 <= m 
                G_right = Givens_rotation([Bk(k, k + 1); Bk(k, k + 2)]);
                Bk(:, k + 1:k + 2) = Bk(:, k + 1:k + 2) * G_right;
                V(:, k + 1:k + 2) = V(:, k + 1:k + 2) * G_right;
            end
            % disp("after right elim superdiagonal")
            % Bk
        end
    end

    Bk
    S = abs(diag(Bk(1:min_dim, 1:min_dim))); %removed abs to match matlab svd
    [S, idx] = sort(S, 'descend');
    U = U(:, 1:length(S));
    V = V(:, 1:length(S));
    U = U(:, idx);
    V = V(:, idx);
end


%Copied from Week 10 Givens_rotation.m

function G = Givens_rotation( x )
    %Givens_rotation Compute Givens rotation G so that G' * x = || x ||_2
    % e_0

    [ m, n ] = size( x );

    assert( m==2 && n==1, 'x must be 2 x 1' );

    normx = norm( x );

    gamma = x(1) / normx;
    sigma = x(2) / normx;

    G = [ (gamma) (-sigma)
          (sigma)  (gamma) ];    
end