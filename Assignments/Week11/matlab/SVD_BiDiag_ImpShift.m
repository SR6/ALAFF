function [S, U, V] = SVD_BiDiag_ImpShift(B)
    % Assumes B is upper bidiagonal
    [n, ~] = size(B);
    max_iter = 2000; % in case it doesn't converge
    eps_mach = eps;

    % Extract diagonals
    d = diag(B);           % Main diagonal
    e = diag(B, 1);        % Superdiagonal

    U = eye(n);
    V = eye(n);

    for k = n:-1:2
        it = 0;
        while abs(e(k-1)) > eps_mach * sqrt(d(k-1)^2 + d(k)^2)
            if it > max_iter
                warning('Maximum iterations reached without full convergence.');
                break;
            end

            % Wilkinson shift
            delta = (d(k-1)^2 - d(k)^2) / 2;
            if delta == 0
                mu = d(k)^2 - abs(e(k-1));
            else
                mu = d(k)^2 - e(k-1)^2 / (delta + sign(delta) * sqrt(delta^2 + e(k-1)^2));
            end

            % Initial Givens rotation to chase the bulge
            x = d(1)^2 - mu;
            y = d(1)*e(1);

            for i = 1:k-1
                % Givens rotation for V (right side)
                [c, s] = givens(x, y);
                G = [c s; -s c];
                V(:, i:i+1) = V(:, i:i+1) * G;

                % Apply G' on B from the right
                temp_d = d(i);
                temp_e = e(i);
                d(i)   = c*temp_d + s*temp_e;
                e(i)   = -s*temp_d + c*temp_e;
                if i < k-1
                    temp = d(i+1);
                    d(i+1) = c*temp;
                    x = s*temp;
                    y = e(i+1);
                end

                % Givens rotation for U (left side)
                [c, s] = givens(d(i), x);
                G = [c s; -s c];
                U(:, i:i+1) = U(:, i:i+1) * G;

                % Apply G on B from the left
                temp_d = d(i);
                temp_e = e(i);
                d(i)   = c*temp_d + s*x;
                e(i)   = c*temp_e + s*y;
                if i < k-1
                    temp = d(i+1);
                    d(i+1) = -s*temp;
                    x = c*temp;
                    y = 0;
                end
            end
            it = it + 1;
        end

        % Zero out the converged off-diagonal
        if abs(e(k-1)) <= eps_mach * sqrt(d(k-1)^2 + d(k)^2)
            e(k-1) = 0;
        end
    end

    % Construct singular values vector
    S = abs(d);

    % Reorder descending and reorder U, V accordingly
    [S, idx] = sort(S, 'descend');
    U = U(:, idx);
    V = V(:, idx);
end

function [c, s] = givens(a, b)
    if b == 0
        c = 1; s = 0;
    else
        if abs(b) > abs(a)
            r = a / b;
            s = 1 / sqrt(1 + r^2);
            c = s * r;
        else
            r = b / a;
            c = 1 / sqrt(1 + r^2);
            s = c * r;
        end
    end
end
