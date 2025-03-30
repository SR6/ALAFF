function y = SparseMvMult(nzA, ir, ic, x)
    % Perform matrix-vector multiplication y = A * x for a sparse matrix A
    
    % Get the size of the vector
    n = length(x);
    
    % Initialize the result vector y
    y = zeros(n, 1);
    
    % Perform the sparse matrix-vector multiplication
    for i = 1:length(nzA)
        row = ir(i);
        col = ic(i);
        value = nzA(i);
        y(row) = y(row) + value * x(col);
    end
end
