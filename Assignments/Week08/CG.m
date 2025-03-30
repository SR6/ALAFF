function [ x, niters ] = CG( A, b, x0 )
% Method to calculate x using conjugated gradient method

tol = 1e-6; %instead of r = 0 which might not happen with flop accumulated error
max_iter = 100000; % prevent infinite loop

niters = 0;
x = x0;
r = b;

while (norm(r) > tol && niters < max_iter)
    if niters == 0
        p = r;
    else
        gamma = (r' * r)/(r_old' * r_old);
        p = r + gamma * p_old;
    end
    
    alpha = (r' * r) / (p' * A*p);
    x_old = x;
    x = x + alpha * p;
    r_old = r;
    r = r - alpha * A * p;
    p_old = p;
    niters = niters + 1;
end
