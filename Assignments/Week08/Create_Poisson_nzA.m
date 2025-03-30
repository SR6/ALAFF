function [nzA,ir,ic] = Create_Poisson_nzA(N)
%Create_Poisson_nzA Creates the compressed row matrix representation of A
assert(N>0, "Bad N input. N must be positive integer");
% Number of grid points
n = (N+2) * (N+2);

% Initialize the output vectors
nzA= []; %non zero array of A
ir = zeros(n+1,1);%row starting index for first non-zero element
ic = []; %columns associated with nzA elements

% Stencil coefficients for 2D Poisson problem (specific here, not generic)
main_diag = 4;
off_diag = -1;

%Counter for non-zero elements, starts at zero because Matlab indexes at 1
nonZeroCount = 0;

for meshPoint = 1:n
    ir(meshPoint) = nonZeroCount+1;
    

    meshPoint;
    row_idx = floor((meshPoint-1)/(N+2)+1); %provides row index start at 1
    col_idx = mod(meshPoint-1,(N+2))+1; %provide col index start at 1

    % Bottom neighbor
    if row_idx > 1 %row - (N+2) > 0
        nzA = [nzA, off_diag];
        ic = [ic, meshPoint - (N+2)];
        nonZeroCount = nonZeroCount + 1;
    end

    % Left neighbor
    if col_idx > 1 % || mod(row, N) ~= 1
        nzA = [nzA, off_diag];
        ic = [ic, meshPoint - 1];
        nonZeroCount = nonZeroCount + 1;
    end

    % Main diagonal
    nzA = [nzA, main_diag];
    ic = [ic, meshPoint];
    nonZeroCount = nonZeroCount + 1;

    % Right neighbor
    if col_idx < (N+2) %|| mod(row, N) ~= 0
        nzA = [nzA, off_diag];
        ic = [ic, meshPoint + 1];
        nonZeroCount = nonZeroCount + 1;
    end

    % Top neighbor
    if row_idx < (N+2) %row + (N+2) <= n
        nzA = [nzA, off_diag];
        ic = [ic, meshPoint + (N+2)];
        nonZeroCount = nonZeroCount + 1;
    end

    % Final entry for ir array
    ir(n + 1) = nonZeroCount + 1;
    
end
