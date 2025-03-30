function [ x, niters ] = Method_of_Steepest_Descent( A, b, x0 )
% Function to implement the method of Steepest Descent
% to solve Ax = b with initial guess x(0)
% Must count the number of iterations prior to meeting tolerance
% for r(k) = 0, which isn't great computing. Prefer to look for
% r(k) <= tolerance 

tol = 1e-6; %instead of r = 0 which might not happen with flop accumulated error
max_iter = 100000; % prevent infinite loop

niters = 0;
x = x0;
r = b - A * x;

while (norm(r) > tol && niters < max_iter)
    p = r;
    q = A*p;
    alpha = (p'*r)/ (p'*q);
    x = x + alpha * p;
    r = r - alpha * q;
    niters = niters + 1;
end

