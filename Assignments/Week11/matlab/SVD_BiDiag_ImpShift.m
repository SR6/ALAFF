function [S, U, V] = SVD_BiDiag_ImpShift(B)
    % Computes the SVD of a bidiagonal matrix B using implicit shift algorithm
    % Input: B - upper bidiagonal matrix
    % Output: S - vector of singular values (sorted in descending order)
    %         U, V - orthogonal matrices such that B = U * diag(S) * V'

    %disp("in impshift NOT test file")
    [n, m] = size(B);
    min_dim = min(n, m);
    U = eye(n);
    V = eye(m);
    Bk = B;
    max_iter = 10000;% * max(n, m);
    defl_tol = 1e-15;
    %defl_tol = 150 * eps(norm(Bk, 'fro'));
    
    for iter = 1:max_iter
        
        %clearing up small values in the off diagonal
        
        for i = 1:min_dim - 1
            if abs(Bk(i, i + 1)) <= defl_tol %tol * (abs(Bk(i, i)) + abs(Bk(i + 1, i + 1)))
                Bk(i, i + 1) = 0;
            end
        end

        % Check for convergence
        off_diag = abs(diag(Bk, 1));
        if all(off_diag < defl_tol) %|| max(off_diag) < tol * norm(Bk, 'fro')
            disp(["reached convergence in ", iter])
            break;
        end

        % Find the unreduced part of the bidiagonal matrix
        q = min_dim;
        while q > 1 && abs(Bk(q - 1, q)) <= defl_tol
            q = q - 1;
        end
        if q <=1
            continue;
        end
        p = 1;
        while p < q - 1 && abs(Bk(p, p + 1)) <= defl_tol
            p = p + 1;
        end
        
        mu = 0;
        if q > 1
            %wilkinsonShift is function at bottom of this file
            %mu = wilkinsonShift(Bk(q-2,q-1),Bk(q-1,q-1),Bk(q-1,q),Bk(q,q));
            mu = wilkinsonShiftAlt(Bk(q-1,q-1),Bk(q-1,q),Bk(q,q));
        end
        
        

        % Initial Givens rotation to create the bulge
        x = Bk(p, p)^2  - mu; %- Bk(q-2,q-1)^2 + Bk(q-1,q-1)^2
        y = Bk(p, p) * Bk(p, p + 1);
        G_right_0 = Givens_rotation([x; y]);
       
        Bk(:, p:p+1) = Bk(:, p:p+1) * G_right_0;
        V(:, p:p+1) = V(:, p:p+1) * G_right_0;

        % Chase the bulge
        %Givens function at the bottom of this file because it was updated
        %from week 10 version, maybe unnecessarily, but it is working and I
        %don't want to change it.
        for k = p:q - 1
            % Left Givens to zero Bk(k+1, k)
            if k < n
                G = Givens_rotation([Bk(k, k); Bk(k+1, k)]);
                Bk(k:min(n,k+1), :) = G' * Bk(k:min(n,k+1), :);
                U(:, k:k+1) = U(:, k:k+1) * G;
            end
            % Right Givens to zero Bk(k, k+2)
            if k < q - 1 && k + 2 <= m
                G = Givens_rotation([Bk(k, k+1); Bk(k, k+2)]);
                Bk(:, k+1:min(m,k+2)) = Bk(:, k+1:min(m,k+2)) * G;
                V(:, k+1:k+2) = V(:, k+1:k+2) * G;
            end
        end

        %cleaning up tiny off diag values
        Bk(abs(Bk) < defl_tol) = 0;
        
    end

    % Bk
    % Bk(abs(Bk) < tol) = 0;
    Bk = diag(diag(Bk));  %Clearing out any -0.000000 stuff
    S_raw = abs(diag(Bk));

    [S, idx] = sort(S_raw, 'descend');
    U = U(:, idx);
    V = V(:, idx);
    %fixing reordering of U relative to signs
    for i = 1:min_dim
        if (U(:,i)' * (B * V(:,i))) < 0
            U(:,i) = -U(:,i);
        end
    end
end


%Copied and updated from Week 10 Givens_rotation.m

function G = Givens_rotation( x )
    %Givens_rotation Compute Givens rotation G so that G' * x = || x ||_2
    % e_0

    normx = hypot(x(1),x(2)); %norm( x ); updated for numerical stability
    %add for NAN issues when normx is very small
    if normx == 0
        disp("in normx == 0")
        G = eye(2);
        return;
    end

    gamma = x(1) / normx;
    sigma = x(2) / normx;

    G = [ (gamma) (-sigma)
          (sigma)  (gamma) ];    
end
function mu = wilkinsonShift(a0,a, b, c)
    x = a0^2 + a^2;
    y = a*b;
    z = b*c;
    delta = (x - z) / 2;

    mu = (z - abs(delta)*(y^2)) / (abs(delta) + sqrt(delta^2 + y^2));
end

function mu = wilkinsonShiftAlt(a, b, c)
    x = a^2 + b^2;
    y = b * c;
    z = c^2;
    delta = (x - z) / 2;

    mu = z - (y^2) / (abs(delta) + sqrt(delta^2 + y^2));
end