function [nzA,ir,ic] = Create_mdim1_Poisson_problem_nzA(N, m)

%Create_Poisson_nzA Creates the compressed row matrix representation of A
% N is related to number of meshpoints (how many intervals
% m is the number of dimensions
assert(N>0, "Bad N input. N must be positive integer");
assert(m>0, "Bad m input. m must be positive integer");

% Number of grid points
n = (N+2)^m;

% Initialize the output vectors
nzA= []; %non zero array of A
ir = zeros(n+1,1);%row starting index for first non-zero element
ic = []; %columns associated with nzA elements

% Stencil coefficients for mDim Poisson problem
main_diag = 2*m;
off_diag = -1;
% number of possible neighbors is 2*m


%Counter for non-zero elements, starts at zero
nonZeroCount = 0;

for meshpoint = 1:n
    ir(meshpoint) = nonZeroCount+1;

     % Process all "backward" neighbors (smaller indices)
    % for dim = 0:m-1
    %     step = (N+2)^dim; % Neighbor spacing in the current dimension
    %     if mod(meshpoint-1, (N+2)^(dim+1)) >= step % Valid backward neighbor
    %         nzA = [nzA, off_diag];
    %         ic = [ic, meshpoint - step];
    %         nonZeroCount = nonZeroCount + 1;
    %     end
    % end
    for dim = m-1:-1:0
        step = (N+2)^dim; % Neighbor spacing in the current dimension
        if mod(meshpoint - 1, (N+2)^(dim+1)) >= step % Valid backward neighbor
            nzA = [nzA, off_diag];
            ic = [ic, meshpoint - step];
            nonZeroCount = nonZeroCount + 1;
        end
    end

    % Insert the main diagonal after all backward neighbors
    nzA = [nzA, main_diag];
    ic = [ic, meshpoint];
    nonZeroCount = nonZeroCount + 1;

    % Process all "forward" neighbors (larger indices)
    for dim = 0:m-1
        step = (N+2)^dim; % Neighbor spacing in the current dimension
        if mod(meshpoint-1, (N+2)^(dim+1)) < (N+2)^(dim+1) - step % Valid forward neighbor
            nzA = [nzA, off_diag];
            ic = [ic, meshpoint + step];
            nonZeroCount = nonZeroCount + 1;
        end
    end

    % Final row pointer entry
    ir(n+1) = nonZeroCount + 1;
end
